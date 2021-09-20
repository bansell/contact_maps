
#adapted from code printed here: plotly_example("shiny", "event_data")

source('~/Dropbox/Jex_Lab/Swapnil_work/contactMaps_2021.R')

library(shiny)
library(plotly)

ui <- fluidPage(
  plotlyOutput("plot", width="800px",height="800px"),
  #verbatimTextOutput("hover"),
  #verbatimTextOutput("click"),
  verbatimTextOutput("brushing"),
  tableOutput('datTable_all')
  #tableOutput('datTable_head'),
  #tableOutput('datTable_tail'),
  #verbatimTextOutput("selecting")
  #verbatimTextOutput("brushed"),
  #verbatimTextOutput("selected")
)

server <- function(input, output, session) {
  
  nms <- row.names(mtcars)
  
  output$plot <- renderPlotly({
    #p <- ggplotly(ggplot(mtcars, aes(x = mpg, y = wt, customdata = nms)) + geom_point())
    p <- ggplotly(myPlot + coord_equal(), tooltip = 'text' )
    
    p %>% 
      layout(dragmode = "select") %>%
      event_register("plotly_selecting")
  })
  
 
  output$brushing <- renderPrint({
    d <- event_data("plotly_brushing") 
    if (is.null(d)) cat("Drag box dimensions (double-click to clear):") else 
      print(c(cat("Drag box dimensions (double-click to clear):\n"), 
      tibble('residue_1_range' = round(d$x,0), 'residue_2_range' = round(d$y,0))))
  })
  
  
  output$datTable_all <- renderTable({
    d <- event_data("plotly_selecting") 
    if (is.null(d)) "Select residues:" else 
      d %>% dplyr::select(-1) %>% 
      left_join(dat_for_plot, by=c('x'='x1','y'='x2') ) %>% distinct() %>% 
      mutate(across(.cols = c(x,y), .fns= ~ as.integer(.))) %>% 
      dplyr::select(-1) %>% arrange(x,y)})
  
  # output$datTable_head <- renderTable({
  #   d <- event_data("plotly_selecting") 
  #   if (is.null(d)) "Selected residues top:" else d %>% head() %>% dplyr::select(-1)})
  # 
  # output$datTable_tail <- renderTable({
  #   d <- event_data("plotly_selecting") 
  #   if (is.null(d)) "Selected residues bottom:" else d %>% tail() %>% dplyr::select(-1)})
  # 
 
}

shinyApp(ui, server) 
