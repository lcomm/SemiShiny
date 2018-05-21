library("rsconnect")
rsconnect::deployApp(appDir = "~/Dropbox/Semicompeting-PS/SemiShiny/", 
                     appFiles = c("app.R", "Strategy1models.md",
                                  "Compartment_shiny.png",
                                  "DESCRIPTION", "README.md"),
                     appName = "Semicompeting-PS",
                     account = "lcomm")

