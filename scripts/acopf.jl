using JuMP, Ipopt;
using PowerSystems;


function run_ACOPF(case_name)


    # Leemos los datos utilizando PowerSystems
    dat = PowerSystems.PowerModelsData(case_name).data

    # Los siguientes únicamente son aliases
    branches = dat["branch"]
    buses = dat["bus"]
    gens = dat["gen"]
    loads = dat["load"]
    shunts = dat["shunt"]

    # Creamos los conjuntos a ser utilizados en los indices de JuMP
    B = keys(buses)
    L = keys(branches)
    G = keys(gens)
    C = keys(loads)
    S = keys(shunts)


    # Mapeo: bus -> [cargas]
    load_map = Dict((i, String[]) for i in B)
    for (key, load) in loads
        push!(load_map[string(load["load_bus"])], key)
    end

    # Mapeo: bus -> [gens]
    gen_map = Dict((i, String[]) for i in keys(buses))
    for (key, load) in gens
        push!(gen_map[string(load["gen_bus"])], key)
    end

    # Mapeo: bus -> [shunts]
    shunt_map = Dict((i, String[]) for i in keys(buses))
    for (key, load) in shunts
        push!(shunt_map[string(load["shunt_bus"])], key)
    end

    # Creamos el modelo de JuMP y asignamos el solver a utilizar (Ipopt)
    mod = Model(Ipopt.Optimizer)


    @variable(mod, buses[k]["vmin"] <= v[k in B] <= buses[k]["vmax"], start = 1.0)

    @variable(mod, δ[k in B], start = 0.0)


    for (k, bus_data) in buses

        # Limites de voltaje
        #set_lower_bound(v[k], bus_data["vmin"]);
        #set_upper_bound(v[k], bus_data["vmax"]);

        # Referencia angular para el nodo "slack"
        if (bus_data["bus_type"] == 3)
            set_lower_bound(δ[k], bus_data["va"])
            set_upper_bound(δ[k], bus_data["va"])
        end
    end

    # Next dictionaries map bus->{Pkm, Qkm, Pmk, Qmk}
    Pinj_f_bra = Dict((k, []) for k in B)
    Pinj_t_bra = Dict((k, []) for k in B)
    Qinj_f_bra = Dict((k, []) for k in B)
    Qinj_t_bra = Dict((k, []) for k in B)


    for brn in values(branches)

        # Nodos k y m
        k = string(brn["f_bus"])
        m = string(brn["t_bus"])

        # Calculamos la admitancia serie entre los nodos k y m
        Ys = 1.0 / (brn["br_r"] + im * brn["br_x"])
        Bc = brn["b_fr"]
        Gc = brn["g_fr"]
        # Calculamos el tap complejo (magnitud de tap y ángulo de desfase)
        cTap = brn["tap"] * exp(deg2rad(brn["shift"]) * im)

        # Admitancia vista desde el nodo m
        Ymm = Ys + Gc + im * Bc
        # Admitancia vista desde el nodo k
        Ykk = Ymm / brn["tap"]^2
        # Admitancia entre los nodos k y m
        Ykm = -Ys / conj(cTap)
        # Admitancia entre los nodos m y k
        Ymk = -Ys / cTap
        brn["Gkk"] = real(Ykk)
        brn["Bkk"] = imag(Ykk)
        brn["Gkm"] = real(Ykm)
        brn["Bkm"] = imag(Ykm)
        brn["Gmm"] = real(Ymm)
        brn["Bmm"] = imag(Ymm)
        brn["Gmk"] = real(Ymk)
        brn["Bmk"] = imag(Ymk)

        if (brn["br_status"] == 0)
            @info(
                "La rama %s entre los nodos %d y %d está inactiva",
                key,
                brn["f_bus"],
                brn["t_bus"]
            )
            continue
        end

        # Expresion no-lineal del flujo de potencia activa y reactiva entre el nodo k y m
        Pkm = @NLexpression(
            mod,
            v[k]^2 * brn["Gkk"] +
            v[k] * v[m] * (brn["Gkm"] * cos(δ[k] - δ[m]) + brn["Bkm"] * sin(δ[k] - δ[m]))
        )
        Qkm = @NLexpression(
            mod,
            -v[k]^2 * brn["Bkk"] +
            v[k] * v[m] * (brn["Gkm"] * sin(δ[k] - δ[m]) - brn["Bkm"] * cos(δ[k] - δ[m]))
        )

        # Expresion no-lineal del flujo de potencia activa y reactiva entre el nodo m y k
        Pmk = @NLexpression(
            mod,
            v[m]^2 * brn["Gmm"] +
            v[m] * v[k] * (brn["Gmk"] * cos(δ[m] - δ[k]) + brn["Bmk"] * sin(δ[m] - δ[k]))
        )
        Qmk = @NLexpression(
            mod,
            -v[m]^2 * brn["Bmm"] +
            v[m] * v[k] * (brn["Gmk"] * sin(δ[m] - δ[k]) - brn["Bmk"] * cos(δ[m] - δ[k]))
        )

        # Incluímos la expresión no-lineal en cada una de las listas de inyecciones a nivel nodo
        push!(Pinj_f_bra[k], Pkm)
        push!(Qinj_f_bra[k], Qkm)
        push!(Pinj_t_bra[m], Pmk)
        push!(Qinj_t_bra[m], Qmk)

        # En caso que tengamos un límite térmico que cuidar. Agregamos la restricción
        if (brn["rate_a"] > 0)
            # Restricción de potencia aparente: Skm <= Smax
            @NLconstraint(mod, Pkm^2 + Qkm^2 <= brn["rate_a"]^2)
            # Restricción de potencia aparnte: Smk <= Smax
            @NLconstraint(mod, Pmk^2 + Qmk^2 <= brn["rate_a"]^2)
        end


    end

    @variable(mod, gens[i]["pmin"] <= Pg[i in G] <= gens[i]["pmax"])
    @variable(mod, gens[i]["qmin"] <= Qg[i in G] <= gens[i]["qmax"])


    # Pg - Pd - Psh - Pkm == 0
    @NLconstraint(mod, dP[k in B], sum(Pg[i] for i in gen_map[k]) -
                                   sum(loads[load]["pd"] for load in load_map[k]) -
                                   sum(shunts[shunt]["gs"] * v[k]^2 for shunt in shunt_map[k]) -
                                   sum(Pkm for Pkm in Pinj_f_bra[k]) -
                                   sum(Pmk for Pmk in Pinj_t_bra[k]) ==0 );
    # Qg - Qd + Qsh - Qkm == 0
    @NLconstraint(mod, dQ[k in B], sum(Qg[i] for i in gen_map[k]) -
                                   sum(loads[load]["qd"] for load in load_map[k]) +
                                   sum(shunts[shunt]["bs"] * v[k]^2 for shunt in shunt_map[k]) -
                                   sum(Qkm for Qkm in Qinj_f_bra[k]) -
                                   sum(Qmk for Qmk in Qinj_t_bra[k]) ==0 );


    @objective(
        mod,
        Min,
        sum(
            sum(gen["cost"][j] * Pg[i]^(gen["ncost"] - j) for j = 1:gen["ncost"])
            for (i, gen) in gens
        )
    )

    optimize!(mod)

    return dat, mod
end
