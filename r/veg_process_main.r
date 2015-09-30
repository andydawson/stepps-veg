run_40knots_c1 = list(suff_data = '12taxa_6341cells_40knots_v0.4',
                      suff_fit  = '12taxa_6341cells_40knots_v0.4_c1')
run_80knots_c1 = list(suff_data = '12taxa_6341cells_80knots_v0.4',
                  suff_fit  = '12taxa_6341cells_80knots_v0.4_c1')
run_120knots_c1 = list(suff_data = '12taxa_6341cells_120knots_v0.4',
                  suff_fit  = '12taxa_6341cells_120knots_v0.4_c1')
run_160knots_c1 = list(suff_data = '12taxa_6341cells_160knots_v0.4',
                  suff_fit  = '12taxa_6341cells_160knots_v0.4_c1')
run_200knots_c1 = list(suff_data = '12taxa_6341cells_200knots_v0.4',
                  suff_fit  = '12taxa_6341cells_200knots_v0.4_c1')
run_220knots_c1 = list(suff_data = '12taxa_6341cells_220knots_v0.4',
                  suff_fit  = '12taxa_6341cells_220knots_v0.4_c1')
run_240knots_c1 = list(suff_data = '12taxa_6341cells_240knots_v0.4',
                       suff_fit  = '12taxa_6341cells_240knots_v0.4_c1')
run_260knots_c1 = list(suff_data = '12taxa_6341cells_260knots_v0.4',
                       suff_fit  = '12taxa_6341cells_260knots_v0.4_c1')


runs = list(run_40knots_c1, 
            run_80knots_c1,
            run_120knots_c1, 
            run_160knots_c1,
            run_200knots_c1, 
            run_220knots_c1,
            run_240knots_c1,
            run_260knots_c1)

# runs = list(run_200knots_c1, run_200knots_c2)

for (run in runs) {
  source('r/veg_process.r')
}