using ArchGDAL
const AG = ArchGDAL
using Infiltrator
# using DataStructures
using DelimitedFiles

function get_vulnerability(data_obj)
    if isfile("data_to_vulnerability.txt")
        # @infiltrate
        vul_mat = readdlm("data_to_vulnerability.txt")
        data_to_vulnerability = Dict(zip(convert(Vector{Int64}, vul_mat[:,1]), vul_mat[:,2]))
    else
        data_to_vulnerability = read_shapefile_and_return_vulnerability(get_coords_from_data(data_obj))
        # @infiltrate
        open("data_to_vulnerability.txt", "w") do io
            writedlm(io, data_to_vulnerability)
        end
    end

    return data_to_vulnerability
end

function get_coords_from_data(data_obj)
    bus_to_coords = Dict()
    for busname in keys(data_obj["bus"])
        # @infiltrate
        bus_data = data_obj["bus"][busname]
        if haskey(bus_data, "lat") && haskey(bus_data, "lon")
            bus_to_coords[busname] = (bus_data["lat"], bus_data["lon"])
        # else
        #     bus_to_coords[busname] = nothing
        end
    end
    return bus_to_coords
end

function read_shapefile_and_return_vulnerability(bus_coord_list, filename="test/data/shapefile")
    data_to_vulnerability = Dict()#DefaultDict(0)
    # AG.registerdrivers() do
        AG.read(filename) do filename_dataset
            # println(AG.getproj(filename_dataset))
            filename_gis = AG.getlayer(filename_dataset, 0)
            # println(AG.getspatialref(filename_gis))
            num_features = AG.nfeature(filename_gis)
            # println(num_features)
            for i in 0:(num_features - 1)
                ArchGDAL.getfeature(filename_gis, i) do feature
                    overall_vulnerability = AG.getfield(feature, "overall_sc")
                    county_name = AG.getfield(feature, "county_nam")
                    geom = AG.getgeom(feature)
                    for (coord_name, coords) in bus_coord_list
                        coords_point = AG.createpoint(coords)
                        source = AG.importEPSG(4326)
                        target = AG.importEPSG(3857)
                        AG.createcoordtrans(source, target) do transform
                            AG.transform!(coords_point, transform)
                            if AG.contains(geom, coords_point)
                                # data_to_vulnerability[coords] = (coord_name, county_name, overall_vulnerability)
                                data_to_vulnerability[coord_name] = overall_vulnerability
                            end
                        end
                    end
                end
            end
        end
    # end
    return data_to_vulnerability
end