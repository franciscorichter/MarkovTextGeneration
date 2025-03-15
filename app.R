library(shiny)
library(bslib)
library(igraph)
library(visNetwork)

#------------------------------------------------------------
# 1) Preprocessing function
#------------------------------------------------------------
preproc <- function(canon) {
  letter  <- "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRTSUVWXYZ"
  letter  <- substring(letter, 1:nchar(letter), 1:nchar(letter))
  delimit <- ".,;:?!"
  delimit <- substring(delimit, 1:nchar(delimit), 1:nchar(delimit))
  remov   <- c("\\", "\"", "/", "'")
  newcanon <- NULL
  n <- length(canon)
  
  for (i in seq_len(n)) {
    nlet <- nchar(canon[i])
    while (nlet > 0 && sum(substr(canon[i], nlet, nlet) == remov) > 0) {
      canon[i] <- substr(canon[i], 1, nlet - 1)
      nlet <- nlet - 1
    }
    if (nchar(canon[i]) > 0 && sum(substr(canon[i], 1, 1) == letter) > 0) {
      nlet <- nchar(canon[i])
      # If last char is punctuation, split it off
      if (sum(substr(canon[i], nlet, nlet) == delimit) > 0) {
        newcanon <- c(newcanon, substr(canon[i], 1, nlet - 1), substr(canon[i], nlet, nlet))
      } else {
        newcanon <- c(newcanon, canon[i])
      }
    }
  }
  newcanon
}

#------------------------------------------------------------
# 2) Markov chain text generation function
#------------------------------------------------------------
writer <- function(canon, n.words = 50, depth = 2, seed = NULL) {
  n.tot <- length(canon)
  delimit <- c(".,;:?!")
  delimit <- substring(delimit, 1:nchar(delimit), 1:nchar(delimit))
  
  if (is.null(seed) || length(seed) == 0) {
    if (n.tot < depth) return("")  # not enough words to seed
    start <- sample(1:(n.tot - depth), 1)
    book <- canon[start:(start + depth - 1)]
  } else {
    book <- seed
  }
  
  while (length(book) < n.words) {
    lb <- length(book)
    tempdepth <- if (lb < depth) lb else depth
    select <- 1:(n.tot - 1)
    for (i in seq_len(tempdepth)) {
      idx <- which(canon == book[lb - i + 1])
      idx <- idx[idx < n.tot]
      idx <- idx + i - 1
      select <- intersect(select, idx)
    }
    if (length(select) == 0) {
      selection <- sample(canon, 1)
    } else if (length(select) == 1) {
      selection <- canon[select + 1]
    } else {
      selection <- sample(canon[select + 1], 1)
    }
    book <- c(book, selection)
  }
  
  verse <- ""
  for (token in book) {
    if (token %in% delimit) {
      verse <- paste0(verse, token)
    } else {
      verse <- paste(verse, token)
    }
  }
  paste(strwrap(verse, width = 80), collapse = "\n")
}

#------------------------------------------------------------
# 3) Transition matrix
#------------------------------------------------------------
transition <- function(canon) {
  words <- unique(canon)
  T.mat <- matrix(0, nrow = length(words), ncol = length(words))
  for (i in seq_along(canon)[-1]) {
    from <- which(words == canon[i - 1])
    to   <- which(words == canon[i])
    T.mat[from, to] <- T.mat[from, to] + 1
  }
  T.mat <- t(apply(T.mat, 1, function(x) if (sum(x) > 0) x / sum(x) else x))
  dimnames(T.mat) <- list(words, words)
  T.mat
}

#------------------------------------------------------------
# 4) Gather subgraph up to `chain_depth` steps from `start_word`
#------------------------------------------------------------
gather_subgraph <- function(T.mat, start_word, chain_depth = 1) {
  all_words <- rownames(T.mat)
  if (! (start_word %in% all_words)) {
    return(list(
      nodes = data.frame(id = "Word not found"),
      edges = data.frame()
    ))
  }
  
  visited <- character(0)    # track visited words
  frontier <- c(start_word)  # current boundary
  edges_df <- data.frame(from=character(0), to=character(0),
                         weight=numeric(0), label=character(0))
  
  for (step in seq_len(chain_depth)) {
    new_frontier <- character(0)
    for (w in frontier) {
      idx <- which(T.mat[w, ] > 0)
      if (length(idx) == 0) next
      successors <- colnames(T.mat)[idx]
      # Create edges for each successor
      for (s in successors) {
        w_prob <- T.mat[w, s]
        edges_df <- rbind(edges_df, data.frame(
          from = w,
          to   = s,
          weight = w_prob,
          label  = sprintf("%.2f", w_prob),
          stringsAsFactors = FALSE
        ))
      }
      new_frontier <- c(new_frontier, successors)
    }
    visited <- unique(c(visited, frontier))
    frontier <- setdiff(unique(new_frontier), visited)
  }
  visited <- unique(c(visited, frontier))
  
  nodes_df <- data.frame(id = visited, label = visited, stringsAsFactors = FALSE)
  edges_df <- edges_df[edges_df$from %in% visited & edges_df$to %in% visited, ]
  
  list(nodes = nodes_df, edges = edges_df)
}

#------------------------------------------------------------
# 5) Mapping display names to file names (in data/ subfolder)
#------------------------------------------------------------
text_files <- c(
  "Shakespeare - All's Well"                 = "allswell.txt",
  "Oscar Wilde - The Picture of Dorian Gray" = "doriangray.txt",
  "Douglas Adams - Hitchhiker's Guide"       = "DouglasAdams-HitchhikersGuide.txt",
  "Jane Austen - Pride and Prejudice"        = "JaneAusten-PrideAndPrejudice.txt",
  "E. E. Cummings - The Enormous Room"       = "EECummings-TheEnormousRoom.txt",
  "Joseph Conrad - Heart of Darkness"        = "JosephConrad-HeartOfDarkness.txt",
  "Aristotle - Categories"                   = "Aristotle-Categories.txt",
  "Lewis Carroll - Alice in Wonderland"      = "LewisCarroll-AlicesAdventuresInWonderland.txt",
  "James Joyce - Ulysses"                    = "JamesJoyce-Ulysses.txt",
  "Bram Stoker - Dracula"                    = "BramStoker-Dracula.txt",
  "Homer - Iliad"                            = "Homer-Iliad.txt",
  "Charles Darwin - On the Origin of Species"= "CharlesDarwin-OnTheOriginOfSpecies.txt",
  "King James Bible"                         = "KingJamesBible.txt",
  "Leo Tolstoy - War and Peace"              = "LeoTolstoy-WarAndPeace.txt",
  "Friedrich Nietzsche - Beyond Good & Evil" = "FriedrichNietzsche-BeyondGoodAndEvil.txt",
  "Charles Dickens - Great Expectations"     = "CharlesDickens-GreatExpectations.txt",
  "Dr Seuss"                                 = "drseuss.txt"
)

#------------------------------------------------------------
# 6) Shiny UI
#------------------------------------------------------------
ui <- fluidPage(
  theme = bs_theme(bootswatch = "minty", version = 5),
  titlePanel("Markov Chain Text Generator & Network"),
  
  tabsetPanel(
    tabPanel("Text Generator",
             sidebarLayout(
               sidebarPanel(
                 selectInput("selected_text", "Select a text source:", names(text_files)),
                 numericInput("n_words", "Number of words to generate:", 50, min = 1, max = 2000),
                 numericInput("depth", "Depth of Markov chain:", 2, min = 1, max = 5),
                 textInput("seed", "Seed text (optional, space-separated):", ""),
                 actionButton("generate", "Generate", class = "btn-primary")
               ),
               mainPanel(
                 verbatimTextOutput("generated_text")
               )
             )
    ),
    tabPanel("Markov Network",
             sidebarLayout(
               sidebarPanel(
                 selectInput("selected_text_network", "Select text for network:", names(text_files)),
                 selectizeInput("starting_word", "Starting word:", choices = NULL,
                                options = list(placeholder = 'Select a starting word')),
                 numericInput("chain_steps", "Chain depth (steps):", 1, min = 1, max = 5),
                 actionButton("show_network", "Show Network", class = "btn-primary")
               ),
               mainPanel(
                 visNetworkOutput("markov_network", height = "600px")
               )
             )
    )
  )
)

#------------------------------------------------------------
# 7) Shiny Server
#------------------------------------------------------------
server <- function(input, output, session) {
  
  # --- Text Generation ---
  observeEvent(input$generate, {
    file_name <- text_files[[input$selected_text]]
    file_path <- file.path("data", file_name)
    
    raw_text <- tryCatch(
      scan(file_path, what = "character", sep = " ", quote = "",
           strip.white = TRUE, blank.lines.skip = TRUE),
      error = function(e) character(0)
    )
    
    if (length(raw_text) == 0) {
      output$generated_text <- renderText("Error: Could not read the file or the file is empty.")
      return()
    }
    
    canon <- preproc(raw_text)
    seed_vec <- if (nchar(input$seed) == 0) NULL else strsplit(input$seed, "\\s+")[[1]]
    result <- writer(canon, n.words = input$n_words, depth = input$depth, seed = seed_vec)
    output$generated_text <- renderText(result)
  })
  
  # --- Update "starting_word" choices on selecting text (Markov Network) ---
  observeEvent(input$selected_text_network, {
    file_name <- text_files[[input$selected_text_network]]
    file_path <- file.path("data", file_name)
    
    raw_text <- tryCatch(
      scan(file_path, what = "character", sep = " ", quote = "",
           strip.white = TRUE, blank.lines.skip = TRUE),
      error = function(e) character(0)
    )
    if (length(raw_text) == 0) return()
    
    canon <- preproc(raw_text)
    words <- sort(unique(canon))
    
    # Populate selectizeInput with all unique words
    updateSelectizeInput(session, "starting_word", choices = words, server = TRUE)
  })
  
  # --- Show Markov Network subgraph from "starting_word" ---
  observeEvent(input$show_network, {
    file_name <- text_files[[input$selected_text_network]]
    file_path <- file.path("data", file_name)
    
    raw_text <- tryCatch(
      scan(file_path, what = "character", sep = " ", quote = "",
           strip.white = TRUE, blank.lines.skip = TRUE),
      error = function(e) character(0)
    )
    
    if (length(raw_text) == 0) {
      output$markov_network <- renderVisNetwork({
        visNetwork(nodes = data.frame(id = "Error"), edges = data.frame())
      })
      return()
    }
    
    canon <- preproc(raw_text)
    T.mat <- transition(canon)
    
    subg <- gather_subgraph(T.mat, start_word = input$starting_word,
                            chain_depth = input$chain_steps)
    nodes <- subg$nodes
    edges <- subg$edges
    
    # If there's an error or no edges
    if (nrow(nodes) == 1 && nodes$id[1] == "Word not found") {
      output$markov_network <- renderVisNetwork({
        visNetwork(nodes = data.frame(id = "Starting word not found"), edges = data.frame())
      })
      return()
    }
    if (nrow(edges) == 0) {
      if (nrow(nodes) == 0) {
        output$markov_network <- renderVisNetwork({
          visNetwork(nodes = data.frame(id = "No words found"), edges = data.frame())
        })
      } else {
        output$markov_network <- renderVisNetwork({
          visNetwork(nodes, edges = data.frame())
        })
      }
      return()
    }
    
    # --- Highlight the starting word in a different color ---
    # Give all nodes a default color
    nodes$color.background <- "#B3CDE0"  # light bluish
    nodes$color.border <- "#00526D"
    
    # Give the starting word a distinct color
    start_idx <- which(nodes$id == input$starting_word)
    if (length(start_idx) == 1) {
      nodes$color.background[start_idx] <- "#FFD700"  # gold
      nodes$color.border[start_idx] <- "#FFA500"      # orange
    }
    
    output$markov_network <- renderVisNetwork({
      visNetwork(nodes, edges) %>%
        visEdges(arrows = "to") %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visLayout(randomSeed = 123)
    })
  })
}

#------------------------------------------------------------
# 8) Run the app
#------------------------------------------------------------
shinyApp(ui, server)
