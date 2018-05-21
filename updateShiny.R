library("rsconnect")
rsconnect::deployApp(appDir = "~/Dropbox/Semicompeting-PS/SemiShiny/", 
                     appFiles = c("app.R", "Strategy1models.md",
                                  "DESCRIPTION", "README.md"),
                     appName = "Semicompeting-PS",
                     account = "lcomm")

