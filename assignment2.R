library(shiny)
library(Seurat)
library(openai)

Sys.setenv(OPENAI_API_KEY = 'XXXXXX')

data("pbmc_small")

# UI for application 
ui <- fluidPage(
  titlePanel("GPTCelltype Annotation"), 
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Select Input Data Format", 
                  choices = list("Input Gene list" = "1", "Differential Gene Table" = "2", "Seurat Object" = "3", "Other" = "4"), 
                  selected = "Input Gene list"),
      conditionalPanel(
        condition = "input.dataset == '1'",
        textAreaInput("gene_list", "Enter Gene List (New line separated for each group. Space separated for each gene)", 
                      "", rows = 5), 
        actionButton("example", "Example Query")
      ),
      conditionalPanel(
        condition = "input.dataset == '2'", 
        fileInput("file", "Upload markers file", accept = c(".rds"))
      ),
      conditionalPanel(
        condition = "input.dataset == '3'", 
        fileInput("file2", "Upload pbmc", accept = c(".rds"))
      ), 
      selectizeInput("model", label = "Select Model", 
                     choices = c("gpt-4o", "gpt-4", "gpt-3.5-turbo"), 
                     selected = "gpt-4o", multiple = FALSE),
      actionButton("annotate", "Annotate Cell Types")
    ), 
    mainPanel(
      titlePanel("Cell Type Annotation: "), 
      conditionalPanel(
        condition = "input.dataset == '1'", 
        tableOutput("table")
      ),
      conditionalPanel(
        condition = "input.dataset == '2'", 
        tableOutput("table3")
      ), 
      conditionalPanel(
        condition = "input.dataset == '3'", 
        tableOutput("table2"), 
        plotOutput("dimPlot")
      ),
      conditionalPanel(
        condition = "input.dataset == '4'", 
        textOutput("other")
      )
    )
  )
)

# Server for application 
server <- function(input, output, session) {
  
  observeEvent(input$annotate, {
    output$dimPlot <- NULL 
    output$other <- NULL
    output$table <- NULL
    
    if (input$dataset == 1) {
      gene_groups <- strsplit(input$gene_list, "\n")[[1]]
      gene_list <- lapply(gene_groups, function(x) strsplit(x, " ")[[1]])
      names(gene_list) <- paste0("Group", seq_along(gene_list))
      
      res <- gptcelltype(input = gene_list, tissuename = NULL, model = input$model, topgenenumber = 10) 
      output$table <- renderTable({
        res_table <- data.frame(Group = names(res), Genes = sapply(gene_list, paste, collapse = " "), Cell_Type = res)
        res_table
      })
      
    } else if (input$dataset == 2) {
      req(input$file)
      
      file_data <- readRDS(input$file$datapath)
      all.markers <- file_data
      
      res <- gptcelltype(all.markers, 
                         tissuename = 'human PBMC', 
                         model = input$model)
      output$table3 <- renderTable({
        res_table3 <- data.frame(Group = names(res), Cell_Type = res)
        res_table3
      })
      
    } else if (input$dataset == 3) {
      req(input$file2)
      
      dataset <- readRDS(input$file2$datapath)
      
      suppressWarnings({
        all.markers <- FindAllMarkers(object = dataset)
      })
      res <- gptcelltype(all.markers, 
                         tissuename = 'human PBMC', 
                         model = input$model)
      dataset@meta.data$celltype <- as.factor(res[as.character(Idents(dataset))])
      
      output$dimPlot <- renderPlot({
        DimPlot(dataset, group.by = 'celltype')
      })
      
      output$table2 <- renderTable({
        res_table2 <- data.frame(Group = names(res), Cell_Type = res)
        res_table2
      })
    } else {
      output$other <- renderText({
        "Sorry we are not supporting this option!!"
      })
    }
  })
  
  observeEvent(input$example, {
    updateTextAreaInput(session, "gene_list", value = "CD4 CD3D\nCD14")
  })
}

# GPT cell type function 
gptcelltype <- function(input, tissuename = NULL, model = 'gpt-4', topgenenumber = 10) {
  OPENAI_API_KEY <- Sys.getenv("OPENAI_API_KEY")
  if (OPENAI_API_KEY == "") {
    print("Note: OpenAI API key not found: returning the prompt itself.")
    API.flag <- 0
  } else {
    API.flag <- 1
  }
  
  if (class(input) == 'list') {
    input <- sapply(input, paste, collapse = ',')
  } else {
    input <- input[input$avg_log2FC > 0,, drop = FALSE]
    input <- tapply(input$gene, list(input$cluster), function(i) paste0(i[1:topgenenumber], collapse = ','))
  }
  
  if (!API.flag) {
    message <- paste0('Identify cell types of ', tissuename, ' cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types. ',  "\n", paste0(names(input), ':', unlist(input), collapse = "\n"))
    return(message)
  } else {
    print("Note: OpenAI API key found: returning the cell type annotations.")
    cutnum <- ceiling(length(input) / 30)
    if (cutnum > 1) {
      cid <- as.numeric(cut(1:length(input), cutnum))	
    } else {
      cid <- rep(1, length(input))
    }
    
    allres <- sapply(1:cutnum, function(i) {
      id <- which(cid == i)
      flag <- 0
      while (flag == 0) {
        k <- openai::create_chat_completion(
          model = model,
          message = list(list("role" = "user", "content" = paste0('Identify cell types of ', tissuename, ' cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types.\n', paste(input[id], collapse = '\n'))))
        )
        res <- strsplit(k$choices[,'message.content'], '\n')[[1]]
        if (length(res) == length(id))
          flag <- 1
      }
      names(res) <- names(input)[id]
      res
    }, simplify = F) 
    print('Note: It is always recommended to check the results returned by GPT-4 in case of\n AI hallucination, before going to down-stream analysis.')
    return(gsub(',$', '', unlist(allres)))
  }
}

# Run the app 
shinyApp(ui, server)