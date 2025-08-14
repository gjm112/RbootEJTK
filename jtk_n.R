jtk_n <- function(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec){
  
  input <- expand_grid(
    phase = phase_vec,
    shift = phase_shift_vec,
    asymmetry = asymmetry_vec
  ) %>%
    filter(asymmetry < phase) 
  
  out <- pmap_dfr(list(input$phase, input$shift, input$asymmetry),
                  ~{
    temp <- jtk1(timepoints, values, ..1, ..2, ..3)
    tibble(
      phase = ..1,
      shift = ..2,
      asymmetry = ..3,
      tau = temp[1],
      pval = temp[2]
    )
    }
  ) %>% arrange(pval)
  
  # results <- data.frame()
  # for (p in phase_vec){
  #   for (s in phase_shift_vec){
  #     for (a in asymmetry_vec){
  #       if (a < p){
  #     temp <- jtk1(timepoints, values, p, s, a)
  #     results <- bind_rows(results, data.frame(phase = p, shift = s , asymmetry = a, tau = temp[1], pval = temp[2]))
  #       }
  #     }
  #   }
  # }
  # results <- results |> arrange(pval)
  return(head(out,1))
}
