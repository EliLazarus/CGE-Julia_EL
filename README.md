# CGE-Julia_EL
An open-source Computable General Equilirium (CGE) model in the Julia language.

This project is inactive and closed in favor of https://github.com/anthofflab/MPSGE.jl . MPSGE.jl is a Julia package to enable simpler and free CGE modelling. The package (will) allow for building models by just structuring the data and providing elasticity parameters.


This is work is now defunct, initially developing toward a simple version of a multi-sector Computable General Equilirium, with taxes.
We were building from basic material from the 2010 EcoMod CGE Modeling workshop, and Cecilia Springer's https://github.com/ccspringer/CGE-in-JuMP/blob/master/Model.jl project.
Following Cecilia, we use the Julia JuMP Package, and the Ipopt solver.
Eventually, this should build to be able to use real data and generate results from experiments. This model will aim to be more or less nuetral in terms of region/nation, and number of sectors etc. The particular goal is to model multiple externality taxes.
