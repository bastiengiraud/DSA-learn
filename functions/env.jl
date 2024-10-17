# install necessary packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.status()
Pkg.add(PackageSpec(name="PowerModels", version="0.21.2"))
Pkg.add(PackageSpec(name="PowerModelsAnnex", version="0.11.0"))
Pkg.add(PackageSpec(name="PowerModelsSecurityConstrained", version="0.12.0"))
Pkg.add(PackageSpec(name="PowerSimulationsDynamics", version="0.14.2"))
Pkg.add(PackageSpec(name="PowerSystemCaseBuilder", version="1.2.5"))
Pkg.add(PackageSpec(name="PowerSystems", version="3.3.0"))
Pkg.add(PackageSpec(name="Sundials", version="4.25.0"))
Pkg.add(PackageSpec(name="Tables", version="1.12.0"))
Pkg.add(PackageSpec(name="TimeSeries", version="0.23.2"))
Pkg.add(PackageSpec(name="LatinHypercubeSampling", version="1.9.0"))
Pkg.add(PackageSpec(name="StatsPlots", version="0.15.7"))
Pkg.add(PackageSpec(name="RCall", version="0.14.6"))
Pkg.add(PackageSpec(name="Plots", version="1.40.8"))
Pkg.add("MosekTools")
Pkg.add("Suppressor")
# ENV["R_HOME"] = "C:\\Program Files\\R\\R-4.3.2"
Pkg.add("Polyhedra")
Pkg.add("LinearAlgebra")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("DelimitedFiles")
Pkg.add("Statistics")
Pkg.add("Distributions")
Pkg.add("Memento")
Pkg.add("InfrastructureModels")
Pkg.add("Ipopt")
Pkg.add("UUIDs")

















# Pkg.update("LatinHypercubeSampling")
# Pkg.update("PowerModelsAnnex")
# Pkg.update("PowerModelsSecurityConstrained")
# Pkg.add("MosekTools")
# Pkg.build("StatsPlots")
# Pkg.add("Suppressor")
# ENV["R_HOME"] = "C:\\Program Files\\R\\R-4.3.2"
# Pkg.build("RCall")
# Pkg.add("Polyhedra")
# Pkg.add("LinearAlgebra")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Statistics")
# Pkg.add("Distributions")
# Pkg.build("Plots")
# Pkg.add("Memento")
# Pkg.add("InfrastructureModels")