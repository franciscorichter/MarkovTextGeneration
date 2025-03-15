library(shiny)
library(bslib)
library(visNetwork)

# -------------------------------------------------------------------------
# 1) Load Preprocessed Data
# -------------------------------------------------------------------------
preprocessed_corpora <- readRDS("preprocessed_corpora.rds")

# -------------------------------------------------------------------------
# 2) Map user-friendly names to corpus keys
# -------------------------------------------------------------------------
text_files <- c(
  "Shakespeare - All's Well"                 = "allswell",
  "Oscar Wilde - The Picture of Dorian Gray" = "doriangray",
  "Douglas Adams - Hitchhiker's Guide"       = "DouglasAdams-HitchhikersGuide",
  "Jane Austen - Pride and Prejudice"        = "JaneAusten-PrideAndPrejudice",
  "E. E. Cummings - The Enormous Room"       = "EECummings-TheEnormousRoom",
  "Joseph Conrad - Heart of Darkness"        = "JosephConrad-HeartOfDarkness",
  "Aristotle - Categories"                   = "Aristotle-Categories",
  "Lewis Carroll - Alice in Wonderland"      = "LewisCarroll-AlicesAdventuresInWonderland",
  "James Joyce - Ulysses"                    = "JamesJoyce-Ulysses",
  "Bram Stoker - Dracula"                    = "BramStoker-Dracula",
  "Homer - Iliad"                            = "Homer-Iliad",
  "Charles Darwin - On the Origin of Species"= "CharlesDarwin-OnTheOriginOfSpecies",
  "King James Bible"                         = "KingJamesBible",
  "Leo Tolstoy - War and Peace"              = "LeoTolstoy-WarAndPeace",
  "Friedrich Nietzsche - Beyond Good & Evil" = "FriedrichNietzsche-BeyondGoodAndEvil",
  "Charles Dickens - Great Expectations"     = "CharlesDickens-GreatExpectations",
  "Dr Seuss"                                 = "drseuss"
)

# -------------------------------------------------------------------------
# 3) Shared Functions
# -------------------------------------------------------------------------

# Simple Markov chain text generator
generate_chain <- function(canon, n.words = 10, depth = 2, seed = NULL) {
  n.tot <- length(canon)
  if (n.tot == 0) return(character(0))
  
  if (is.null(seed) || length(seed) == 0) {
    if (n.tot < depth) return(character(0))
    start <- sample(1:(n.tot - depth), 1)
    book  <- canon[start:(start + depth - 1)]
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
  book
}

# Build transition matrix
build_transition_matrix <- function(canon) {
  words <- unique(canon)
  T.mat <- matrix(0, nrow = length(words), ncol = length(words))
  for (i in seq_along(canon)[-1]) {
    from <- which(words == canon[i - 1])
    to   <- which(words == canon[i])
    T.mat[from, to] <- T.mat[from, to] + 1
  }
  T.mat <- t(apply(T.mat, 1, function(x) if (sum(x) > 0) x/sum(x) else x))
  dimnames(T.mat) <- list(words, words)
  T.mat
}

# Build a chain tree with top-K alternatives at each step.
build_chain_tree <- function(chain, T.mat, topK = 3) {
  if (length(chain) < 1) {
    nodes <- data.frame(id = 1, label = "Empty chain", stringsAsFactors = FALSE)
    edges <- data.frame()
    return(list(nodes = nodes, edges = edges))
  }
  
  words <- rownames(T.mat)
  nodes <- data.frame(id = integer(0), label = character(0), stringsAsFactors = FALSE)
  nodes$color.background <- character(0)
  nodes$color.border <- character(0)
  
  edges <- data.frame(from = integer(0), to = integer(0),
                      label = character(0), color = character(0),
                      stringsAsFactors = FALSE)
  
  node_counter <- 0
  create_node <- function(word, step, is_chosen = FALSE) {
    label_text <- paste0(step, ": ", word)
    node_counter <<- node_counter + 1
    new_id <- node_counter
    if (is_chosen) {
      bg_col <- "#FFD700"  # gold
      br_col <- "#FFA500"  # orange
    } else {
      bg_col <- "#D1E1FF"
      br_col <- "#666666"
    }
    nodes <<- rbind(nodes, data.frame(
      id = new_id, label = label_text,
      color.background = bg_col, color.border = br_col,
      stringsAsFactors = FALSE
    ))
    return(new_id)
  }
  
  create_edge <- function(from_id, to_id, prob, chosen = FALSE) {
    edge_col <- if (chosen) "red" else "#999999"
    lbl <- sprintf("%.2f", prob)
    edges <<- rbind(edges, data.frame(
      from = from_id, to = to_id,
      label = lbl, color = edge_col, stringsAsFactors = FALSE
    ))
  }
  
  current_id <- create_node(chain[1], 1, is_chosen = TRUE)
  
  for (i in seq_len(length(chain) - 1)) {
    current_word <- chain[i]
    next_word <- chain[i + 1]
    step_number <- i + 1
    if (!(current_word %in% words)) break
    row_idx <- which(words == current_word)
    row_probs <- T.mat[row_idx, ]
    nonzero_idx <- which(row_probs > 0)
    if (length(nonzero_idx) == 0) break
    sorted_idx <- nonzero_idx[order(row_probs[nonzero_idx], decreasing = TRUE)]
    keep_idx <- head(sorted_idx, topK)
    
    chosen_child_id <- NA
    for (idx in keep_idx) {
      candidate_word <- words[idx]
      prob_val <- row_probs[idx]
      is_chosen_child <- (candidate_word == next_word)
      child_id <- create_node(candidate_word, step_number, is_chosen = is_chosen_child)
      create_edge(current_id, child_id, prob_val, chosen = is_chosen_child)
      if (is_chosen_child) chosen_child_id <- child_id
    }
    if (is.na(chosen_child_id)) break
    current_id <- chosen_child_id
  }
  list(nodes = nodes, edges = edges)
}

# Gather a subgraph for the Markov network
gather_subgraph <- function(T.mat, start_word, chain_steps = 1) {
  all_words <- rownames(T.mat)
  if (!(start_word %in% all_words)) {
    return(list(
      nodes = data.frame(id = "Word not found", label = "Word not found", stringsAsFactors = FALSE),
      edges = data.frame()
    ))
  }
  visited <- character(0)
  frontier <- c(start_word)
  edges_df <- data.frame(from = character(0), to = character(0),
                         weight = numeric(0), label = character(0),
                         stringsAsFactors = FALSE)
  for (step in seq_len(chain_steps)) {
    new_frontier <- character(0)
    for (w in frontier) {
      idx <- which(T.mat[w, ] > 0)
      if (length(idx) == 0) next
      successors <- colnames(T.mat)[idx]
      for (s in successors) {
        w_prob <- T.mat[w, s]
        edges_df <- rbind(edges_df, data.frame(
          from = w, to = s, weight = w_prob,
          label = sprintf("%.2f", w_prob),
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

# -------------------------------------------------------------------------
# 4) UI: Single Sidebar for All Tabs
# -------------------------------------------------------------------------
ui <- fluidPage(
  theme = bs_theme(
    bootswatch = "cosmo",
    base_font = font_google("Roboto"),
    code_font = font_google("Fira Code")
  ),
  titlePanel("Markov Text Generator"),
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_text", "Select a text source:", names(text_files)),
      numericInput("n_words", "Number of words to generate:", 10, min = 1, max = 2000),
      numericInput("depth", "Markov depth:", 2, min = 1, max = 5),
      selectizeInput("seed_words", "Seed words (optional):", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "Type or select seed words", create = TRUE)),
      # Extra parameters for Tree and Network tabs:
      conditionalPanel(
        condition = "input.main_tab != 'Text Generator'",
        hr(),
        h4("Additional Parameters")
      ),
      conditionalPanel(
        condition = "input.main_tab == 'Chain Tree'",
        numericInput("top_k_tree", "Top K alternatives:", 3, min = 1, max = 10)
      ),
      conditionalPanel(
        condition = "input.main_tab == 'Markov Network'",
        selectizeInput("starting_word", "Starting word:", choices = NULL,
                       options = list(placeholder = "Select or type a word", create = TRUE)),
        numericInput("chain_steps", "Network chain depth:", 1, min = 1, max = 10)
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "main_tab",
        tabPanel("Text Generator",
                 br(),
                 actionButton("generate_text", "Generate Text", class = "btn-primary"),
                 br(), br(),
                 verbatimTextOutput("generated_text")
        ),
        tabPanel("Chain Tree",
                 br(),
                 actionButton("generate_tree", "Generate Chain Tree", class = "btn-primary"),
                 br(), br(),
                 visNetworkOutput("chain_tree_network", height = "600px")
        ),
        tabPanel("Markov Network",
                 br(),
                 actionButton("generate_network", "Generate Markov Network", class = "btn-primary"),
                 br(), br(),
                 visNetworkOutput("markov_network", height = "600px")
        )
      )
    )
  )
)

# -------------------------------------------------------------------------
# 5) Server
# -------------------------------------------------------------------------
server <- function(input, output, session) {
  # Update seed words for Text Generator when text is selected
  observeEvent(input$selected_text, {
    corpus_key <- text_files[[input$selected_text]]
    canon <- preprocessed_corpora[[corpus_key]]
    if (is.null(canon) || length(canon) == 0) {
      updateSelectizeInput(session, "seed_words", choices = NULL, server = TRUE)
    } else {
      all_words <- sort(unique(canon))
      updateSelectizeInput(session, "seed_words", choices = all_words, server = TRUE)
    }
  })
  
  # -----------------------------
  # Tab 1: Text Generator
  # -----------------------------
  observeEvent(input$generate_text, {
    corpus_key <- text_files[[input$selected_text]]
    canon <- preprocessed_corpora[[corpus_key]]
    if (is.null(canon) || length(canon) == 0) {
      output$generated_text <- renderText("Error: This corpus is empty or not found.")
      return()
    }
    seed_vec <- input$seed_words
    if (length(seed_vec) == 0) seed_vec <- NULL
    chain <- generate_chain(canon, n.words = input$n_words, depth = input$depth, seed = seed_vec)
    if (length(chain) == 0) {
      output$generated_text <- renderText("No text generated (not enough words?).")
    } else {
      output$generated_text <- renderText({
        paste(chain, collapse = " ")
      })
    }
  })
  
  # -----------------------------
  # Tab 2: Chain Tree
  # -----------------------------
  observeEvent(input$generate_tree, {
    corpus_key <- text_files[[input$selected_text]]
    canon <- preprocessed_corpora[[corpus_key]]
    if (is.null(canon) || length(canon) == 0) {
      output$chain_tree_network <- renderVisNetwork({
        visNetwork(nodes = data.frame(id = "Error"), edges = data.frame())
      })
      return()
    }
    T.mat <- build_transition_matrix(canon)
    seed_vec <- input$seed_words
    if (length(seed_vec) == 0) seed_vec <- NULL
    chain <- generate_chain(canon, n.words = input$n_words, depth = input$depth, seed = seed_vec)
    if (length(chain) < 1) {
      output$chain_tree_network <- renderVisNetwork({
        visNetwork(nodes = data.frame(id = "No words generated"), edges = data.frame())
      })
      return()
    }
    subg <- build_chain_tree(chain, T.mat, topK = input$top_k_tree)
    nodes <- subg$nodes
    edges <- subg$edges
    
    # Safeguard: if the tree is too big, ask user to lower parameters
    if (nrow(nodes) > 100) {
      output$chain_tree_network <- renderVisNetwork({
        visNetwork(nodes = data.frame(id = "Too many nodes; please lower chain depth or top K."), edges = data.frame())
      })
      return()
    }
    
    output$chain_tree_network <- renderVisNetwork({
      visNetwork(nodes, edges) %>%
        visEdges(arrows = "to", smooth = FALSE) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visHierarchicalLayout(direction = "LR", levelSeparation = 200)
    })
  })
  
  # -----------------------------
  # Tab 3: Markov Network
  # -----------------------------
  observeEvent(input$selected_text, {  # update starting word choices when text changes
    corpus_key <- text_files[[input$selected_text]]
    canon <- preprocessed_corpora[[corpus_key]]
    if (is.null(canon) || length(canon) == 0) {
      updateSelectizeInput(session, "starting_word", choices = NULL, server = TRUE)
    } else {
      all_words <- sort(unique(canon))
      updateSelectizeInput(session, "starting_word", choices = all_words, server = TRUE)
    }
  })
  
  observeEvent(input$generate_network, {
    corpus_key <- text_files[[input$selected_text]]
    canon <- preprocessed_corpora[[corpus_key]]
    if (is.null(canon) || length(canon) == 0) {
      output$markov_network <- renderVisNetwork({
        visNetwork(nodes = data.frame(id = "Error"), edges = data.frame())
      })
      return()
    }
    T.mat <- build_transition_matrix(canon)
    subg <- gather_subgraph(T.mat, input$starting_word, input$chain_steps)
    nodes <- subg$nodes
    edges <- subg$edges
    if (nrow(nodes) == 1 && nodes$id[1] == "Word not found") {
      output$markov_network <- renderVisNetwork({
        visNetwork(nodes = data.frame(id = "Starting word not found", label = "Starting word not found"), edges = data.frame())
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
    nodes$color.background <- "#D1E1FF"
    nodes$color.border <- "#666666"
    start_idx <- which(nodes$id == input$starting_word)
    if (length(start_idx) == 1) {
      nodes$color.background[start_idx] <- "#FFD700"
      nodes$color.border[start_idx] <- "#FFA500"
    }
    output$markov_network <- renderVisNetwork({
      visNetwork(nodes, edges) %>%
        visEdges(arrows = "to") %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visLayout(randomSeed = 123)
    })
  })
}

shinyApp(ui, server)

