refVals = []
refPrecs = []

#! format: off
# SC ========== Test number 1 reference values and precision match template. =======
# SC ========== /home/jmc/cliMa/cliMa_new_jmc/test/Ocean/SplitExplicit/simple_dbl_gyre.jl test reference values ======================================
# BEGIN SCPRINT
# varr - reference values (from reference run)
# parr - digits match precision (hand edit as needed)
#
# [
#  [ MPIStateArray Name, Field Name, Maximum, Minimum, Mean, Standard Deviation ],
#  [         :                :          :        :      :          :           ],
# ]
varr = [
 [  "oce Q_3D",   "u[1]", -3.14619475384711300148e-01,  3.39267297882365259643e-01,  7.01853683660589589000e-02,  6.26660169588590110212e-02 ],
 [  "oce Q_3D",   "u[2]", -5.40115064220653584570e-01,  5.60182007126385217433e-01, -1.33501073697162020437e-02,  6.39836694061109295673e-02 ],
 [  "oce Q_3D",       :η, -4.06675526509356544125e+00,  3.31166018225353031568e+00, -3.26112220789978460647e-04,  2.19172771408384203795e+00 ],
 [  "oce Q_3D",       :θ, -2.74364012693550339550e-02,  2.48197520207026265382e+01,  8.97230380533195415182e+00,  6.53046685376568358805e+00 ],
 [   "oce aux",       :w, -2.03640396387641661388e-03,  1.27826397939490195299e-03, -1.32450141281447979874e-06,  1.89396215557024257155e-04 ],
 [   "oce aux",    :pkin, -6.47709914591248292481e+01,  0.00000000000000000000e+00, -1.29101328536586201778e+01,  1.55140959870780044128e+01 ],
 [   "oce aux",     :wz0, -2.20202997677884930623e-04,  1.75117332904405536903e-04, -7.41190918600093764232e-10,  1.05290510005729925458e-04 ],
 [   "oce aux", "u_d[1]", -3.38043452790656284712e-01,  2.93237463017143928923e-01,  3.32762033529279827038e-02,  5.96482633726020602949e-02 ],
 [   "oce aux", "u_d[2]", -4.61068840293685944243e-01,  6.19043054353931254674e-01,  7.51573507616540206072e-03,  6.25045068699146971758e-02 ],
 [   "oce aux", "ΔGu[1]", -4.59854094097925271542e-05,  6.52335864749383005954e-06, -9.46878317787071804603e-07,  4.75766999399048653958e-06 ],
 [   "oce aux", "ΔGu[2]", -9.58414599150476549379e-06,  1.47783639396590521823e-05,  6.55418429718145221163e-06,  3.22275048104942662288e-06 ],
 [   "oce aux",       :y,  0.00000000000000000000e+00,  6.00000000000000093132e+06,  3.00000000000000000000e+06,  1.73273876310492726043e+06 ],
 [ "baro Q_2D",   "U[1]", -2.68001600158847033128e+01,  3.54678958344688567195e+02,  1.10718760115483163986e+02,  7.02151015509662528302e+01 ],
 [ "baro Q_2D",   "U[2]", -4.58800625138165628414e+02,  3.08927975166181340683e+02, -6.24392179443026691388e+01,  6.58632718135061310250e+01 ],
 [ "baro Q_2D",       :η, -4.06651472360696875086e+00,  3.30995829157504939388e+00, -3.26054200343120720394e-04,  2.19300873645531391176e+00 ],
 [  "baro aux",  "Gᵁ[1]", -1.95700759424814914322e-02,  1.37956228229377586558e-01,  2.84063495336121622781e-03,  1.42734794289876652101e-02 ],
 [  "baro aux",  "Gᵁ[2]", -4.43350918189771570077e-02,  2.87524379745142943943e-02, -1.96625528915443520406e-02,  9.66856943716619572637e-03 ],
 [  "baro aux", "U_c[1]", -2.66901779856414016479e+01,  3.54311741174972667068e+02,  1.10718668449473625515e+02,  7.01384607022182819946e+01 ],
 [  "baro aux", "U_c[2]", -4.59179604966002727906e+02,  3.08121811474068351799e+02, -6.26075213657350246876e+01,  6.58906196442101048660e+01 ],
 [  "baro aux",     :η_c, -4.06675526509356544125e+00,  3.31166018225353031568e+00, -3.26112220790006758324e-04,  2.19179980119349959722e+00 ],
 [  "baro aux", "U_s[1]", -2.66900837860213968611e+01,  3.54313639706073615798e+02,  1.10719880610984176883e+02,  7.01390166119913089915e+01 ],
 [  "baro aux", "U_s[2]", -4.59180095304110864163e+02,  3.08121734125270336335e+02, -6.26079379020178237170e+01,  6.58908018110790578703e+01 ],
 [  "baro aux",     :η_s, -4.06676620312277403713e+00,  3.31166762340400744336e+00, -3.26113359381027595116e-04,  2.19180135714922696977e+00 ],
 [  "baro aux",  "Δu[1]", -3.25140056764470369594e-04,  3.73797414826839259974e-04,  1.99435503229933815453e-05,  1.05093079175937200950e-04 ],
 [  "baro aux",  "Δu[2]", -9.06671179261403877752e-05,  5.96317424935956285140e-04,  1.05449620472310954004e-04,  1.37755025647565924191e-04 ],
 [  "baro aux",  :η_diag, -4.06374752927915761092e+00,  3.30838227635390991210e+00, -3.23174912520092548832e-04,  2.19135433800063950116e+00 ],
 [  "baro aux",      :Δη, -3.41173917778858637462e-03,  4.32125937995797571034e-03, -2.93730826990848609067e-06,  1.07007365822442868271e-03 ],
 [  "baro aux",       :y,  0.00000000000000000000e+00,  6.00000000000000093132e+06,  3.00000000000000000000e+06,  1.73279575381979602389e+06 ],
]
parr = [
 [  "oce Q_3D",   "u[1]",    12,    12,    12,    12 ],
 [  "oce Q_3D",   "u[2]",    12,    12,    12,    12 ],
 [  "oce Q_3D",       :η,    12,    12,     8,    12 ],
 [  "oce Q_3D",       :θ,    12,    12,    12,    12 ],
 [   "oce aux",       :w,    12,    12,     8,    12 ],
 [   "oce aux",    :pkin,    12,    12,    12,    12 ],
 [   "oce aux",     :wz0,    12,    12,     8,    12 ],
 [   "oce aux", "u_d[1]",    12,    12,    12,    12 ],
 [   "oce aux", "u_d[2]",    12,    12,    12,    12 ],
 [   "oce aux", "ΔGu[1]",    12,    12,    12,    12 ],
 [   "oce aux", "ΔGu[2]",    12,    12,    12,    12 ],
 [   "oce aux",       :y,    12,    12,    12,    12 ],
 [ "baro Q_2D",   "U[1]",    12,    12,    12,    12 ],
 [ "baro Q_2D",   "U[2]",    12,    12,    12,    12 ],
 [ "baro Q_2D",       :η,    12,    12,     8,    12 ],
 [  "baro aux",  "Gᵁ[1]",    12,    12,    12,    12 ],
 [  "baro aux",  "Gᵁ[2]",    12,    12,    12,    12 ],
 [  "baro aux", "U_c[1]",    12,    12,    12,    12 ],
 [  "baro aux", "U_c[2]",    12,    12,    12,    12 ],
 [  "baro aux",     :η_c,    12,    12,     8,    12 ],
 [  "baro aux", "U_s[1]",    12,    12,    12,    12 ],
 [  "baro aux", "U_s[2]",    12,    12,    12,    12 ],
 [  "baro aux",     :η_s,    12,    12,     8,    12 ],
 [  "baro aux",  "Δu[1]",    12,    12,    12,    12 ],
 [  "baro aux",  "Δu[2]",    12,    12,    12,    12 ],
 [  "baro aux",  :η_diag,    12,    12,     8,    12 ],
 [  "baro aux",      :Δη,    12,    12,     8,    12 ],
 [  "baro aux",       :y,    12,    12,    12,    12 ],
]
# END SCPRINT
# SC ====================================================================================

    append!(refVals ,[ varr ] )
    append!(refPrecs,[ parr ] )

#! format: on