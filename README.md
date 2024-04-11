# R studio

These analyses are part of my paper entitled: Marquez, V., Carbone, L.M., Chiapero, A.L., Ashworth, L., Calviño, A., Zamudio, F., & Aguilar, R. 2022. Silvopastoral and peasant management effects on vegetation and soil quality in the arid Chaco of central Argentina. Journal of Arid Environments. 206,104845. https://doi.org/10.1016/j.jaridenv.2022.104845

To compare the similarity in species composition between different conditions of land use type (silvopastoral and peasant) we built a matrix with values of Bray-Curtis distance from abundance species data of each sampling plot differentiated by strata (herb, shrub and tree). Based on this, we conducted a one-way non-parametric similarity analysis, ANOSIM (with 999 permutations), to determine whether species composition differed between land management. Also, we performed NMDS (Non-Metric Multidimensional Scaling) ordination analysis using the calculated dissimilarity measures to visually assess differences in species composition between management types. This graphical output shows similarity in species composition among parcels, whereby the closer they are ordinated the more similar their composition. We calculated the stress value that represents the difference between distance in the reduced dimension ordination and the complete multidimensional space, thus, it reflects how well the ordination summarizes the observed distances among the parcels. According to Clarke (1993), a stress value lower than 0.1 gives a “good ordination with no real risk of drawing false inferences”. These analyses were run using the vegan package (Oksanen et al., 2016)
 
