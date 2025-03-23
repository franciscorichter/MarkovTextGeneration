library(shiny)
library(bslib)
library(visNetwork)

# ----------------------------------------------------------------------------
# 1) Load Preprocessed Data
# ----------------------------------------------------------------------------
preprocessed_corpora <- readRDS("preprocessed_corpora.rds")

# ----------------------------------------------------------------------------
# 2) Map User-Friendly Names to Corpus Keys
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# 3) Markov Chain Generator for Text
# ----------------------------------------------------------------------------
generate_chain <- function(canon, n.words = NULL, n.phrases = 2,depth = 2, seed = NULL) {
  n.tot <- length(canon)
  if (n.tot == 0) return(character(0))
  
  # If no seed, pick a random start
  if (is.null(seed) || length(seed) == 0) {
    if (n.tot < depth) return(character(0))
    start <- sample(1:(n.tot - depth), 1)
    while(canon[start] %in% c(".", ",", ";", ":", "?", "!")){
      start <- sample(1:(n.tot - depth), 1)
    }
    book  <- canon[start:(start + depth - 1)]
  } else {
    book <- seed
  }
  
  
  if(is.null(n.words)){
    # Build chain
    while (length(which(book %in% c(".", "?", "!")))<n.phrases) {
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

  }else{
  
  # Build chain
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
    
  }
  
  book
}

# ----------------------------------------------------------------------------
# 4) Build Transition Matrix (for the tree)
# ----------------------------------------------------------------------------
build_transition_matrix <- function(canon) {
  words <- unique(canon)
  T.mat <- matrix(0, nrow = length(words), ncol = length(words))
  for (i in seq_along(canon)[-1]) {
    from <- which(words == canon[i - 1])
    to   <- which(words == canon[i])
    T.mat[from, to] <- T.mat[from, to] + 1
  }
  # Convert counts to probabilities
  T.mat <- t(apply(T.mat, 1, function(x) {
    if (sum(x) > 0) x / sum(x) else x
  }))
  dimnames(T.mat) <- list(words, words)
  T.mat
}

# ----------------------------------------------------------------------------
# 5) Build Chain Tree
#     - n_evolutions = how many words in the chain
#     - branching_factor = K, the max # of alt edges per node
#     - We do NOT break if the actual chosen next word is outside top K.
# ----------------------------------------------------------------------------
build_chain_tree <- function(chain, T.mat, branching_factor = 3) {
  
  # If chain is empty, return an "empty" network
  if (length(chain) < 1) {
    nodes <- data.frame(id = 1, label = "Empty chain", stringsAsFactors = FALSE)
    edges <- data.frame()
    return(list(nodes = nodes, edges = edges))
  }
  
  words <- rownames(T.mat)
  
  # Data frames for nodes & edges
  nodes <- data.frame(
    id                = integer(0),
    label            = character(0),
    color.background = character(0),
    color.border     = character(0),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from  = integer(0),
    to    = integer(0),
    label = character(0),
    color = character(0),
    stringsAsFactors = FALSE
  )
  
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
      id                = new_id,
      label            = label_text,
      color.background = bg_col,
      color.border     = br_col,
      stringsAsFactors = FALSE
    ))
    return(new_id)
  }
  
  create_edge <- function(from_id, to_id, prob, chosen = FALSE) {
    edge_col <- if (chosen) "red" else "#999999"
    lbl      <- sprintf("%.2f", prob)
    edges <<- rbind(edges, data.frame(
      from  = from_id,
      to    = to_id,
      label = lbl,
      color = edge_col,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create a node for the first word in the chain
  current_id <- create_node(chain[1], step = 1, is_chosen = TRUE)
  
  # For each step i in [1, length(chain)-1]
  for (i in seq_len(length(chain) - 1)) {
    current_word <- chain[i]
    next_word    <- chain[i + 1]
    step_number  <- i + 1
    
    if (!(current_word %in% words)) {
      # If current_word doesn't exist in T.mat, we can't proceed
      break
    }
    
    row_idx   <- which(words == current_word)
    row_probs <- T.mat[row_idx, ]
    nonzero_idx <- which(row_probs > 0)
    if (length(nonzero_idx) == 0) {
      break
    }
    
    # Sort them descending by probability
    sorted_idx <- nonzero_idx[order(row_probs[nonzero_idx], decreasing = TRUE)]
    # Keep top K
    keep_idx   <- head(sorted_idx, branching_factor)
    
    # Ensure the chosen next_word is also displayed, even if it's outside top K
    chosen_idx <- which(words == next_word)
    if (length(chosen_idx) == 1 && !(chosen_idx %in% keep_idx)) {
      keep_idx <- c(keep_idx, chosen_idx)
    }
    keep_idx <- unique(keep_idx)
    
    chosen_child_id <- NA
    
    # Create child nodes for each candidate in keep_idx
    for (idx in keep_idx) {
      candidate_word   <- words[idx]
      prob_val         <- row_probs[idx]
      is_chosen_child  <- (candidate_word == next_word)
      
      child_id <- create_node(
        word      = candidate_word,
        step      = step_number,
        is_chosen = is_chosen_child
      )
      create_edge(
        from_id = current_id,
        to_id   = child_id,
        prob    = prob_val,
        chosen  = is_chosen_child
      )
      
      if (is_chosen_child) {
        chosen_child_id <- child_id
      }
    }
    
    # Even if the chosen word wasn't in keep_idx, we do NOT break.
    # We continue with the chain. If chosen_child_id is NA, we'll just
    # stay on the same node or break, depending on your preference.
    if (!is.na(chosen_child_id)) {
      current_id <- chosen_child_id
    } else {
      # If there's truly no child for the chosen word, you might break or pick another child.
      # We'll break here, but you could remove this if you want to keep going anyway.
      break
    }
  }
  
  list(nodes = nodes, edges = edges)
}

# ----------------------------------------------------------------------------
# 6) Shiny UI
# ----------------------------------------------------------------------------
ui <- fluidPage(
  theme = bs_theme(
    bootswatch = "cosmo",
    base_font  = font_google("Roboto"),
    code_font  = font_google("Fira Code")
  ),
  
  titlePanel("Markov Text Generator !!"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_text", "Select a text source:", names(text_files)),
      
      # Shared seed words
      selectizeInput(
        "seed_words",
        "Seed words (optional):",
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Type or select seed words",
          create = TRUE
        )
      ),
      
      # Text Generator inputs
      conditionalPanel(
        condition = "input.main_tab == 'Text Generator'",
        
        selectInput("word_or_phrase", "Generate by number of", c("sentences", "words")),
        
        conditionalPanel(
          condition = "input.word_or_phrase == 'words' ",
          numericInput("n_words", "Number of words (Text Gen):", 10, min = 1, max = 2000)
        )
        ,
        
        conditionalPanel(
          condition = "input.word_or_phrase == 'sentences' ",
          numericInput("n_phrases", "Number of sentences (Text Gen):", 2, min = 1, max = 20)
        ),
        
        
        numericInput("markov_depth", "Markov depth (Text Gen):", 2, min = 1, max = 5),
        actionButton("generate_text", "Generate Text", class = "btn-primary")
      ),
      
      # Chain Tree inputs
      conditionalPanel(
        condition = "input.main_tab == 'Chain Tree'",
        numericInput("n_evolutions", "Number of evolutions (Tree length):", 10, min = 1, max = 200),
        numericInput("branching_factor", "Branching factor (Tree width):", 3, min = 1, max = 10),
        actionButton("generate_tree", "Generate Tree", class = "btn-primary")
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "main_tab",
        
        tabPanel("Text Generator",
                 br(),
                 verbatimTextOutput("generated_text")
        ),
        
        tabPanel("Chain Tree",
                 br(),
                 visNetworkOutput("chain_tree_network", height = "600px")
        )
      )
    )
  )
)

# ----------------------------------------------------------------------------
# 7) Shiny Server
# ----------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Update seed words after user picks a text
  observeEvent(input$selected_text, {
    corpus_key <- text_files[[ input$selected_text ]]
    canon <- preprocessed_corpora[[ corpus_key ]]
    
    if (is.null(canon) || length(canon) == 0) {
      updateSelectizeInput(session, "seed_words", choices = NULL, server = TRUE)
    } else {
      all_words <- sort(setdiff(unique(canon), c(".", ",", ";", ":", "?", "!")))
      updateSelectizeInput(session, "seed_words", choices = all_words, server = TRUE)
    }
  })
  
  
  # ---------------------- Tab 1: Text Generator ----------------------
  observeEvent(input$generate_text, {
    corpus_key <- text_files[[ input$selected_text ]]
    canon <- preprocessed_corpora[[ corpus_key ]]
    if (is.null(canon) || length(canon) == 0) {
      output$generated_text <- renderText("Error: This corpus is empty or not found.")
      return()
    }
    
    
    
    seed_vec <- input$seed_words
    if (length(seed_vec) == 0) seed_vec <- NULL
    
    if(input$word_or_phrase=="sentences"){
      n_phrases<-input$n_phrases
      n_words <- NULL
    }else if (input$word_or_phrase=="words"){
      n_words<-input$n_words
      n_phrases <- NULL
    }
    
    
    chain <- generate_chain(
      canon   = canon,
      n.words = n_words,
      n.phrases = n_phrases,
      depth   = input$markov_depth,
      seed    = seed_vec
    )
    
    if (length(chain) == 0) {
      output$generated_text <- renderText("No text generated (maybe not enough words?).")
    } else {
      output$generated_text <- renderText({
        paste(chain, collapse = " ")
      })
    }
  })
  
  # ---------------------- Tab 2: Chain Tree ----------------------
  observeEvent(input$generate_tree, {
    corpus_key <- text_files[[ input$selected_text ]]
    canon <- preprocessed_corpora[[ corpus_key ]]
    if (is.null(canon) || length(canon) == 0) {
      output$chain_tree_network <- renderVisNetwork({
        visNetwork(nodes = data.frame(id="Error"), edges=data.frame())
      })
      return()
    }
    
    # We'll fix Markov depth = 1 for the chain used in the tree
    T.mat <- build_transition_matrix(canon)
    
    seed_vec <- input$seed_words
    if (length(seed_vec) == 0) seed_vec <- NULL
    
    # Generate a chain of length n_evolutions, ignoring top K at generation time
    chain <- generate_chain(
      canon   = canon,
      n.words = input$n_evolutions,
      depth   = 1,   # fixed for the tree
      seed    = seed_vec
    )
    
    if (length(chain) < 1) {
      output$chain_tree_network <- renderVisNetwork({
        visNetwork(nodes = data.frame(id="No words generated"), edges=data.frame())
      })
      return()
    }
    
    subg <- build_chain_tree(
      chain,
      T.mat,
      branching_factor = input$branching_factor
    )
    nodes <- subg$nodes
    edges <- subg$edges
    
    # If the tree is huge, show a message
    if (nrow(nodes) > 150) {
      output$chain_tree_network <- renderVisNetwork({
        visNetwork(
          nodes = data.frame(id="Tree is too large. Lower 'Number of evolutions' or 'Branching factor'."),
          edges = data.frame()
        )
      })
      return()
    }
    
    output$chain_tree_network <- renderVisNetwork({
      visNetwork(nodes, edges) %>%
        visEdges(arrows="to", smooth=FALSE) %>%
        visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE) %>%
        visHierarchicalLayout(direction="LR", levelSeparation=200)
    })
  })
}

shinyApp(ui, server)

