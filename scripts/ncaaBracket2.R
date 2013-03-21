library( randomForest )
library( ggplot2 )
library( grid )

stopQuietly <- function(...)
  {
  blankMsg <- sprintf( "\r%s\r", paste( rep(" ", getOption( "width" ) - 1L ), collapse = " ") );
  stop( simpleError( blankMsg ) );
  } # stopQuietly()

args <- commandArgs( trailingOnly = TRUE )

# inputDir <- args[1]
inputDir <- "~/Desktop"

featureFile <- paste0( inputDir, "/ncaaFeatureTable.RData" )

featureTable <- data.frame()

if( ! is.na( file.info( featureFile )$size ) )
  {
  load( featureFile )
  } else {

  ##################################################################
  #
  # Read data tables from ESPN
  #   1. Read postseason into variable 'featureTable' and keep "team",
  #      "season", and "total points" (we're going to regress on total
  #      points assuming that is a good indicator of who wins the
  #      tournament.)
  #   2. Read regular season stats and fill in featureTable for the
  #      teams that were read for the postseason.  It only makes
  #      sense to keep the regular season stats for teams/seasons
  #      for which postseason stats are available since we're trying
  #      to predict postseason stats based on regular season stats
  #   3. Read AP and USAToday polls
  #
  ##################################################################

  library( "RCurl" )
  library( "XML" )

  cat( "The file ncaaFeatureTable.RData does not exist.  Loading from espn.com\n" )

  baseUrl <- "http://espn.go.com/mens-college-basketball/statistics/team/_/stat/"
  years <- 2002:2012

  ## Read postseason

  pages <- 1:3
  postSeasonUrl <- paste0( baseUrl, "scoring/sort/points/year/" )

  cat( "Reading postseason\n" );
  for( y in years )
    {
    cat( "  ", y, ": page ", sep = "" )
    for( p in pages )
      {
      cat( p )
      count <- ( p - 1 ) * 40 + 1

      pageTable <- readHTMLTable( paste0( postSeasonUrl, y, "/seasontype/3/count/", count ) )
      pageTable <- pageTable$`NULL`
      pageTable <- pageTable[which( pageTable$V2 != 'PER GAME' & pageTable$V2 != 'TEAM'),]

      pageTable <- transform( pageTable, V2 = as.character( V2 ) )
      pageTable <- transform( pageTable, V4 = as.numeric( as.character( V4 ) ) )

      numberOfTeams <- length( pageTable$V2 )

      localFeatureTable <- data.frame( Team=rep.int( 0, numberOfTeams ),
                                       Season=rep.int( 0, numberOfTeams ),
                                       PostSeasonTotalPoints=rep.int( 0, numberOfTeams ) )

      localFeatureTable$Team <- pageTable$V2
      localFeatureTable$Season <- y
      localFeatureTable$PostSeasonTotalPoints <- pageTable$V4

      if( y == years[1] && p == pages[1] )
        {
        featureTable <- localFeatureTable
        }
      else
        {
        featureTable <- rbind( featureTable, localFeatureTable )
        }

      if( p == pages[length( pages )] )
        {
        cat( "\n" )
        }
      else
        {
        cat( ", " )
        }
      }
    }

  ## Read regular seasons

  years <- 2002:2012

  pages <- 1:3

  tableTypes <- c( "avgPoints", "rebounds", "field-goals", "free-throws", "3-points", "assists", "steals", "blocks" )
  headers <- list()
  headers[[1]] <- c( "RK", "TEAM", "GP", "PTS", "FGM-FGA", "FGPercentage_1", "3PM-3PA", "3PPercentage_1", "FTM-FTA", "FTPercentage_1" )
  headers[[2]] <- c( "RK", "TEAM",	"GP",	"OFF",	"ORPG",	"DEF",	"DRPG",	"REB",	"RPG" )
  headers[[3]] <- c( "RK", "TEAM",	"GP",	"PPG_3",	"PerGameFGM",	"PerGameFGA",	"TotalFGM",	"TotalFGA",	"FGPercentage_3",	"2PM_3",	"2PA_3",	"2PPercentage_3",	"PPS_3",	"AdjustedFGPercentage" )
  headers[[4]] <- c( "RK", "TEAM",	"GP",	"PPG_4",	"PerGameFTM",	"PerGameFTA",	"TotalFTM",	"TotalFTA",	"FTPercentage_4" )
  headers[[5]] <- c( "RK",	"TEAM",	"GP",	"PPG_5",	"PerGame3PM",	"PerGame3PA",	"Total3PM",	"Total3PA",	"3PPercentage_5",	"2PM_5",	"2PA_5",	"2PPercentage_5",	"PPS_5",	"FGPercentage_5" )
  headers[[6]] <- c( "RK",	"TEAM",	"GP",	"AST",	"APG",	"TO_6",	"TOPG_6",	"AST/TO" )
  headers[[7]] <- c( "RK",	"TEAM",	"GP",	"STL",	"STPG",	"TO_7",	"TOPG_7",	"PF_7", "ST/TO", "ST/PF" )
  headers[[8]] <- c( "RK",	"TEAM",	"GP",	"BLK",	"PF_8",	"BLKPG",	"BLK/PF" )

  urls <- list()
  for( i in 1:length( tableTypes ) )
    {
    urls[[i]] <- paste0( baseUrl, tableTypes[i], "/sort/" )
    }

  columnsToKeep <- c( "GamesPlayed",                   # = headers[[1]]$GP
                      "PointsPerGame",                 # = headers[[1]]$PTS
                      "FieldGoalPercentage",           # = headers[[1]]$FGPercentage_1
                      "FreeThrowPercentage",           # = headers[[1]]$FTPercentage_1
                      "TotalOffensiveRebounds",        # = headers[[2]]$OFF
                      "OffensiveReboundsPerGame",      # = headers[[2]]$ORPG
                      "TotalDefensiveRebounds",        # = headers[[2]]$DEF
                      "DefensiveReboundsPerGame",      # = headers[[2]]$DRPG
                      "TotalRebounds",                 # = headers[[2]]$REB
                      "ReboundsPerGame",               # = headers[[2]]$RPG
                      "FieldGoalsMadePerGame",         # = headers[[3]]$PerGameFGM
                      "FieldGoalsAttemptedPerGame",    # = headers[[3]]$PerGameFGA
                      "TotalFieldGoalsMade",           # = headers[[3]]$TotalFGM
                      "TotalFieldGoalsAttempted",      # = headers[[3]]$TotalFGA
                      "TwoPointersMade",               # = headers[[3]]$2PM_3
                      "TwoPointersAttempted",          # = headers[[3]]$2PA_3
                      "TwoPointersPercentage",         # = headers[[3]]$2PPercentage_3
                      "PointsPerShot",                 # = headers[[3]]$PPS_3
                      "AdjustedFieldGoalPercentage",   # = headers[[3]]$AdjustedFGPercentage
                      "FreeThrowsMadePerGame",         # = headers[[4]]$PerGameFGM
                      "FreeThrowsAttemptedPerGame",    # = headers[[4]]$PerGameFGA
                      "TotalFreeThrowsMade",           # = headers[[4]]$TotalFTM
                      "TotalFreeThrowsAttempted",      # = headers[[4]]$TotalFTA
                      "ThreePointersMadePerGame",      # = headers[[5]]$PerGame3PM
                      "ThreePointersAttemptedPerGame", # = headers[[5]]$PerGame3PA
                      "TotalThreePointersMade",        # = headers[[5]]$Total3PM
                      "TotalThreePointersAttempted",   # = headers[[5]]$Total3PA
                      "ThreePointPercentage",          # = headers[[5]]$3PPercentage_5
                      "TotalAssists",                  # = headers[[6]]$AST
                      "AssistsPerGame",                # = headers[[6]]$APG
                      "TotalTurnovers",                # = headers[[6]]$TO_6
                      "TurnoversPerGame",              # = headers[[6]]$TOPG_6
                      "AssistsToTurnoversRatio",             # = headers[[6]]$AST/TO
                      "TotalSteals",                   # = headers[[7]]$STL
                      "StealsPerGame",                 # = headers[[7]]$STPG
                      "PersonalFouls",                 # = headers[[7]]$PF_7
                      "StealsToTurnoversRatio",              # = headers[[7]]$ST/TO
                      "StealsToPersonalFoulsRatio",          # = headers[[7]]$ST/PF
                      "Blocks",                        # = headers[[8]]$BLK
                      "BlocksPerGame",                 # = headers[[8]]$BLKPG
                      "BlocksToPersonalFouls"           # = headers[[8]]$BLK/PF
                    )
  featureTable[columnsToKeep] <- 0

  columnsToKeepIndices <- matrix( 0, nrow = length( columnsToKeep ), ncol = 2 )

  columnsToKeepIndices[1,] <- c( 1, 3 )
  columnsToKeepIndices[2,] <- c( 1, 4 )
  columnsToKeepIndices[3,] <- c( 1, 6 )
  columnsToKeepIndices[4,] <- c( 1, 10 )
  columnsToKeepIndices[5,] <- c( 2, 4 )
  columnsToKeepIndices[6,] <- c( 2, 5 )
  columnsToKeepIndices[7,] <- c( 2, 6 )
  columnsToKeepIndices[8,] <- c( 2, 7 )
  columnsToKeepIndices[9,] <- c( 2, 8 )
  columnsToKeepIndices[10,] <- c( 2, 9 )

  columnsToKeepIndices[11,] <- c( 3, 5 )
  columnsToKeepIndices[12,] <- c( 3, 6 )
  columnsToKeepIndices[13,] <- c( 3, 7 )
  columnsToKeepIndices[14,] <- c( 3, 8 )
  columnsToKeepIndices[15,] <- c( 3, 10 )
  columnsToKeepIndices[16,] <- c( 3, 11 )
  columnsToKeepIndices[17,] <- c( 3, 12 )
  columnsToKeepIndices[18,] <- c( 3, 13 )
  columnsToKeepIndices[19,] <- c( 3, 14 )

  columnsToKeepIndices[20,] <- c( 4, 5 )
  columnsToKeepIndices[21,] <- c( 4, 6 )
  columnsToKeepIndices[22,] <- c( 4, 7 )
  columnsToKeepIndices[23,] <- c( 4, 8 )

  columnsToKeepIndices[24,] <- c( 5, 5 )
  columnsToKeepIndices[25,] <- c( 5, 6 )
  columnsToKeepIndices[26,] <- c( 5, 7 )
  columnsToKeepIndices[27,] <- c( 5, 8 )
  columnsToKeepIndices[28,] <- c( 5, 9 )

  columnsToKeepIndices[29,] <- c( 6, 4 )
  columnsToKeepIndices[30,] <- c( 6, 5 )
  columnsToKeepIndices[31,] <- c( 6, 6 )
  columnsToKeepIndices[32,] <- c( 6, 7 )
  columnsToKeepIndices[33,] <- c( 6, 8 )

  columnsToKeepIndices[34,] <- c( 7, 4 )
  columnsToKeepIndices[35,] <- c( 7, 5 )
  columnsToKeepIndices[36,] <- c( 7, 8 )
  columnsToKeepIndices[37,] <- c( 7, 9 )
  columnsToKeepIndices[38,] <- c( 7, 10 )

  columnsToKeepIndices[39,] <- c( 8, 4 )
  columnsToKeepIndices[40,] <- c( 8, 6 )
  columnsToKeepIndices[41,] <- c( 8, 7 )

  for( h in 1:length( headers ) )
    {
    cat( "Reading ", tableTypes[h], "\n", sep = "" )
    for( y in years )
      {
      cat( "  regular season ", y, ": page ", sep = "" )

      for( p in pages )
        {
        cat( p )
        count <- ( p - 1 ) * 40 + 1
        pageTable <- readHTMLTable( paste0( urls[[h]], "/year/", y, "/count/", count ) )
        pageTable <- pageTable$`NULL`
        pageTable <- pageTable[which( pageTable$V2 != 'PER GAME' & pageTable$V2 != 'TEAM'),]

        colnames( pageTable ) <- headers[[h]]
        pageTable <- data.frame( lapply( pageTable, as.character ), stringsAsFactors = FALSE )

        for( j in 1:nrow( pageTable ) )
          {
          idx <- which( featureTable$Team == pageTable$TEAM[j] & featureTable$Season == y )
          if( length( idx ) > 0 )
            {
            for( k in 1:nrow( columnsToKeepIndices ) )
              {
              if( columnsToKeepIndices[k,1] == h )
                {
                featureTable[idx,k+3] <- as.numeric( as.character( pageTable[j,columnsToKeepIndices[k,2]] ) )
                }
              }
            }
          }

        if( p == pages[length( pages )] )
          {
          cat( "\n" )
          }
        else
          {
          cat( ", " )
          }
        }
      }
    }

  ##################################################################
  #
  # Read poll tables from ESPN
  #
  ##################################################################

  years <- 2002:2012
  weeks <- 1:18

  apPolls <- list()
  usaTodayPolls <- list()

  length( apPolls ) <- nrow( featureTable )
  length( usaTodayPolls ) <- nrow( featureTable )

  for( i in 1:length( apPolls ) )
    {
    apPolls[[i]] <- rep( 0, length( weeks ) )
    usaTodayPolls[[i]] <- rep( 0, length( weeks ) )
    }

  cat( "Reading USA Today and AP poll data\n" )

  for( y in years )
    {
    cat( "  Season ", y, ": week ", sep = "" )
    for( w in weeks )
      {
      cat( w )
      pollUrl <- paste0( "http://espn.go.com/mens-college-basketball/rankings/_/", "year/", y, "/week/", w, "/seasontype/2" )
      pageTables <- readHTMLTable( pollUrl )

      # AP Poll

      pageTables[[1]] <- subset( pageTables[[1]], select = c( "V2", "V4" ) )
      pageTables[[1]] <- pageTables[[1]][3:length( pageTables[[1]][,1]),]
      colnames( pageTables[[1]] ) <- c( "Team", "POINTS" )

      pageTables[[1]] <- data.frame( lapply( pageTables[[1]], as.character ), stringsAsFactors = FALSE )

      for( i in 1:length( pageTables[[1]]$Team ) )
        {
        q <- strsplit( as.character( pageTables[[1]]$Team[i] ), "\\(" )[[1]]
        if( length( q ) > 1 )
          {
          r <- strsplit( q[2], "\\)" )[[1]][1]
          if( !is.na( as.numeric( r ) ) )
            {
            pageTables[[1]]$Team[i] <- gsub( "\\s", "", q[1] )
            }
          }
        idx <- which( featureTable$Team == pageTables[[1]]$Team[i] & featureTable$Season == y )
        if( length( idx ) > 0 )
          {
          apPolls[[idx]][w] <- as.numeric( gsub( ",", "", pageTables[[1]]$POINTS[i] ) )
          }
        }

      # USA Today poll

      pageTables[[2]] <- subset( pageTables[[2]], select = c( "V2", "V4" ) )
      pageTables[[2]] <- pageTables[[2]][3:length( pageTables[[2]][,1]),]
      colnames( pageTables[[2]] ) <- c( "Team", "POINTS" )

      pageTables[[2]] <- data.frame( lapply( pageTables[[2]], as.character ), stringsAsFactors = FALSE )

      for( i in 1:length( pageTables[[2]]$Team ) )
        {
        q <- strsplit( as.character( pageTables[[2]]$Team[i] ), "\\(" )[[1]]
        if( length( q ) > 1 )
          {
          r <- strsplit( q[2], "\\)" )[[1]][1]
          if( !is.na( as.numeric( r ) ) )
            {
            pageTables[[2]]$Team[i] <- gsub( "\\s", "", q[1] )
            }
          }
        idx <- which( featureTable$Team == pageTables[[2]]$Team[i] & featureTable$Season == y )
        if( length( idx ) > 0 )
          {
          usaTodayPolls[[idx]][w] <- as.numeric( gsub( ",", "", pageTables[[2]]$POINTS[i] ) )
          }
        }

      if( w == weeks[length(weeks)] )
        {
        cat( "\n", sep = "" )
        }
      else
        {
        cat( ", ", sep = "" )
        }
      }
    }

  pollHeaders <- c( "AP.AverageVotes.1to9weeks", "AP.AverageVotes.10to18weeks", "AP.Trend.13to18weeks",
                    "USAToday.AverageVotes.1to9weeks", "USAToday.AverageVotes.10to18weeks", "USAToday.Trend.13to18weeks"
                  )
  featureTable[pollHeaders] <- 0
  for( i in 1:nrow( featureTable ) )
    {
    featureTable$AP.AverageVotes.1to9weeks[i] <- mean( apPolls[[i]][1:9], na.rm = TRUE )
    featureTable$AP.AverageVotes.10to18weeks[i] <- mean( apPolls[[i]][10:18], na.rm = TRUE )
    featureTable$USAToday.AverageVotes.1to9weeks[i] <- mean( usaTodayPolls[[i]][1:9], na.rm = TRUE )
    featureTable$USAToday.AverageVotes.10to18weeks[i] <- mean( usaTodayPolls[[i]][10:18], na.rm = TRUE )

    spl <- smooth.spline( weeks, apPolls[[i]], spar = 0.5 )
    spl.deriv <- predict( spl, weeks, deriv = 1 )
    featureTable$AP.Trend.13to18weeks[i] <- mean( spl.deriv$y[13:18], na.rm = TRUE )

    spl <- smooth.spline( weeks, usaTodayPolls[[i]], spar = 0.5 )
    spl.deriv <- predict( spl, weeks, deriv = 1 )
    featureTable$USAToday.Trend.13to18weeks[i] <- mean( spl.deriv$y[13:18], na.rm = TRUE )
    }

  # debugging code to look at the spline fits
  #
  # library( ggplot2 )
  #
  # splinePoints <- apPolls[[3]]
  # n <- length( splinePoints )
  # x <- 1:n
  # xx <- seq( 1, n, length.out = 200 )
  # y <- splinePoints
  # d <- data.frame( x = x, y = y )
  # spl <- smooth.spline( x, y, spar = 0.5 )
  # spline.data <- data.frame( y = predict( spl, xx ) )
  # spline.data.deriv <- data.frame( y = predict( spl, xx, deriv = 1 ) )
  # ggplot( d, aes( x, y ) ) + geom_point() +
  #   geom_line( aes( y.x, y.y ), spline.data ) + geom_line( aes( y.x, y.y ), spline.data.deriv )

  ##################################################################
  #
  # Save the feature table for future use
  #
  ##################################################################

  save( featureTable, file = featureFile )
  }

##################################################################
#
# load 2013 regular season
#
##################################################################

teams2013 <- c( 'Louisville', 'North Carolina A&T', 'Colorado State', 'Missouri',
                      'Oklahoma State', 'Oregon', 'Saint Louis', 'New Mexico State',
                      'Memphis', 'Saint Mary\'s', 'Michigan State', 'Valparaiso',
                      'Creighton', 'Cincinnati', 'Duke', 'Albany', ### Midwest
                      'Kansas', 'Western Kentucky', 'North Carolina', 'Villanova',
                      'Virginia Commonwealth', 'Akron', 'Michigan', 'South Dakota State', 'UCLA',
                      'Minnesota', 'Florida', 'Northwestern State',
                      'San Diego State', 'Oklahoma', 'Georgetown',
                      'Florida Gulf Coast',  ### South
                      'Indiana',
                      'James Madison', 'North Carolina State', 'Temple',
                      'UNLV', 'California', 'Syracuse', 'Montana', 'Butler',
                      'Bucknell', 'Marquette', 'Davidson', 'Illinois',
                      'Colorado', 'Miami (FL)', 'Pacific', ### East
                      'Gonzaga', 'Southern University', 'Pittsburgh', 'Wichita State',
                      'Wisconsin', 'Mississippi State', 'Kansas State', 'Boise State', 'La Salle',
                      'Arizona', 'Belmont', 'New Mexico', 'Harvard', 'Notre Dame', 'Iowa State',
                      'Ohio State', 'Iona' ### West
                      )

library( "RCurl" )
library( "XML" )

featureTable2013 <- data.frame( Team = teams2013 )

baseUrl <- "http://espn.go.com/mens-college-basketball/statistics/team/_/stat/"
years <- 2013

pages <- 1:3

## Read regular seasons

pages <- 1:9

tableTypes <- c( "scoring", "rebounds", "field-goals", "free-throws", "3-points", "assists", "steals", "blocks" )
headers <- list()
headers[[1]] <- c( "RK", "TEAM", "GP", "PTS", "FGM-FGA", "FGPercentage_1", "3PM-3PA", "3PPercentage_1", "FTM-FTA", "FTPercentage_1" )
headers[[2]] <- c( "RK", "TEAM",	"GP",	"OFF",	"ORPG",	"DEF",	"DRPG",	"REB",	"RPG" )
headers[[3]] <- c( "RK", "TEAM",	"GP",	"PPG_3",	"PerGameFGM",	"PerGameFGA",	"TotalFGM",	"TotalFGA",	"FGPercentage_3",	"2PM_3",	"2PA_3",	"2PPercentage_3",	"PPS_3",	"AdjustedFGPercentage" )
headers[[4]] <- c( "RK", "TEAM",	"GP",	"PPG_4",	"PerGameFTM",	"PerGameFTA",	"TotalFTM",	"TotalFTA",	"FTPercentage_4" )
headers[[5]] <- c( "RK",	"TEAM",	"GP",	"PPG_5",	"PerGame3PM",	"PerGame3PA",	"Total3PM",	"Total3PA",	"3PPercentage_5",	"2PM_5",	"2PA_5",	"2PPercentage_5",	"PPS_5",	"FGPercentage_5" )
headers[[6]] <- c( "RK",	"TEAM",	"GP",	"AST",	"APG",	"TO_6",	"TOPG_6",	"AST/TO" )
headers[[7]] <- c( "RK",	"TEAM",	"GP",	"STL",	"STPG",	"TO_7",	"TOPG_7",	"PF_7", "ST/TO", "ST/PF" )
headers[[8]] <- c( "RK",	"TEAM",	"GP",	"BLK",	"PF_8",	"BLKPG",	"BLK/PF" )

urls <- list()
for( i in 1:length( tableTypes ) )
  {
  urls[[i]] <- paste0( baseUrl, tableTypes[i], "/seasontype/" )
  }

columnsToKeep <- c( "GamesPlayed",                   # = headers[[1]]$GP
                    "PointsPerGame",                 # = headers[[1]]$PTS
                    "FieldGoalPercentage",           # = headers[[1]]$FGPercentage_1
                    "FreeThrowPercentage",           # = headers[[1]]$FTPercentage_1
                    "TotalOffensiveRebounds",        # = headers[[2]]$OFF
                    "OffensiveReboundsPerGame",      # = headers[[2]]$ORPG
                    "TotalDefensiveRebounds",        # = headers[[2]]$DEF
                    "DefensiveReboundsPerGame",      # = headers[[2]]$DRPG
                    "TotalRebounds",                 # = headers[[2]]$REB
                    "ReboundsPerGame",               # = headers[[2]]$RPG
                    "FieldGoalsMadePerGame",         # = headers[[3]]$PerGameFGM
                    "FieldGoalsAttemptedPerGame",    # = headers[[3]]$PerGameFGA
                    "TotalFieldGoalsMade",           # = headers[[3]]$TotalFGM
                    "TotalFieldGoalsAttempted",      # = headers[[3]]$TotalFGA
                    "TwoPointersMade",               # = headers[[3]]$2PM_3
                    "TwoPointersAttempted",          # = headers[[3]]$2PA_3
                    "TwoPointersPercentage",         # = headers[[3]]$2PPercentage_3
                    "PointsPerShot",                 # = headers[[3]]$PPS_3
                    "AdjustedFieldGoalPercentage",   # = headers[[3]]$AdjustedFGPercentage
                    "FreeThrowsMadePerGame",         # = headers[[4]]$PerGameFGM
                    "FreeThrowsAttemptedPerGame",    # = headers[[4]]$PerGameFGA
                    "TotalFreeThrowsMade",           # = headers[[4]]$TotalFTM
                    "TotalFreeThrowsAttempted",      # = headers[[4]]$TotalFTA
                    "ThreePointersMadePerGame",      # = headers[[5]]$PerGame3PM
                    "ThreePointersAttemptedPerGame", # = headers[[5]]$PerGame3PA
                    "TotalThreePointersMade",        # = headers[[5]]$Total3PM
                    "TotalThreePointersAttempted",   # = headers[[5]]$Total3PA
                    "ThreePointPercentage",          # = headers[[5]]$3PPercentage_5
                    "TotalAssists",                  # = headers[[6]]$AST
                    "AssistsPerGame",                # = headers[[6]]$APG
                    "TotalTurnovers",                # = headers[[6]]$TO_6
                    "TurnoversPerGame",              # = headers[[6]]$TOPG_6
                    "AssistsToTurnoversRatio",             # = headers[[6]]$AST/TO
                    "TotalSteals",                   # = headers[[7]]$STL
                    "StealsPerGame",                 # = headers[[7]]$STPG
                    "PersonalFouls",                 # = headers[[7]]$PF_7
                    "StealsToTurnoversRatio",              # = headers[[7]]$ST/TO
                    "StealsToPersonalFoulsRatio",          # = headers[[7]]$ST/PF
                    "Blocks",                        # = headers[[8]]$BLK
                    "BlocksPerGame",                 # = headers[[8]]$BLKPG
                    "BlocksToPersonalFouls"           # = headers[[8]]$BLK/PF
                  )
featureTable2013[columnsToKeep] <- 0

columnsToKeepIndices <- matrix( 0, nrow = length( columnsToKeep ), ncol = 2 )

columnsToKeepIndices[1,] <- c( 1, 3 )
columnsToKeepIndices[2,] <- c( 1, 4 )
columnsToKeepIndices[3,] <- c( 1, 6 )
columnsToKeepIndices[4,] <- c( 1, 10 )
columnsToKeepIndices[5,] <- c( 2, 4 )
columnsToKeepIndices[6,] <- c( 2, 5 )
columnsToKeepIndices[7,] <- c( 2, 6 )
columnsToKeepIndices[8,] <- c( 2, 7 )
columnsToKeepIndices[9,] <- c( 2, 8 )
columnsToKeepIndices[10,] <- c( 2, 9 )

columnsToKeepIndices[11,] <- c( 3, 5 )
columnsToKeepIndices[12,] <- c( 3, 6 )
columnsToKeepIndices[13,] <- c( 3, 7 )
columnsToKeepIndices[14,] <- c( 3, 8 )
columnsToKeepIndices[15,] <- c( 3, 10 )
columnsToKeepIndices[16,] <- c( 3, 11 )
columnsToKeepIndices[17,] <- c( 3, 12 )
columnsToKeepIndices[18,] <- c( 3, 13 )
columnsToKeepIndices[19,] <- c( 3, 14 )

columnsToKeepIndices[20,] <- c( 4, 5 )
columnsToKeepIndices[21,] <- c( 4, 6 )
columnsToKeepIndices[22,] <- c( 4, 7 )
columnsToKeepIndices[23,] <- c( 4, 8 )

columnsToKeepIndices[24,] <- c( 5, 5 )
columnsToKeepIndices[25,] <- c( 5, 6 )
columnsToKeepIndices[26,] <- c( 5, 7 )
columnsToKeepIndices[27,] <- c( 5, 8 )
columnsToKeepIndices[28,] <- c( 5, 9 )

columnsToKeepIndices[29,] <- c( 6, 4 )
columnsToKeepIndices[30,] <- c( 6, 5 )
columnsToKeepIndices[31,] <- c( 6, 6 )
columnsToKeepIndices[32,] <- c( 6, 7 )
columnsToKeepIndices[33,] <- c( 6, 8 )

columnsToKeepIndices[34,] <- c( 7, 4 )
columnsToKeepIndices[35,] <- c( 7, 5 )
columnsToKeepIndices[36,] <- c( 7, 8 )
columnsToKeepIndices[37,] <- c( 7, 9 )
columnsToKeepIndices[38,] <- c( 7, 10 )

columnsToKeepIndices[39,] <- c( 8, 4 )
columnsToKeepIndices[40,] <- c( 8, 6 )
columnsToKeepIndices[41,] <- c( 8, 7 )

for( h in 1:length( headers ) )
  {
  cat( "Reading ", tableTypes[h], "\n", sep = "" )
  for( y in years )
    {
    cat( "  regular season ", y, ": page ", sep = "" )

    for( p in pages )
      {
      cat( p )
      count <- ( p - 1 ) * 40 + 1

      pageTable <- readHTMLTable( paste0( urls[[h]], "2/count/", count ) )
      pageTable <- pageTable$`NULL`
      pageTable <- pageTable[which( pageTable$V2 != 'PER GAME' & pageTable$V2 != 'TEAM'),]

      colnames( pageTable ) <- headers[[h]]
      pageTable <- data.frame( lapply( pageTable, as.character ), stringsAsFactors = FALSE )

      for( j in 1:nrow( pageTable ) )
        {
        idx <- which( featureTable2013$Team == pageTable$TEAM[j] )
        if( length( idx ) > 0 )
          {
          for( k in 1:nrow( columnsToKeepIndices ) )
            {
            if( columnsToKeepIndices[k,1] == h )
              {
              featureTable2013[idx,k+1] <- as.numeric( as.character( pageTable[j,columnsToKeepIndices[k,2]] ) )
              }
            }
          }
        }

      if( p == pages[length( pages )] )
        {
        cat( "\n" )
        }
      else
        {
        cat( ", " )
        }
      }
    }
  }

##################################################################
#
# Read poll tables from ESPN
#
##################################################################

years <- 2013
weeks <- 1:18

apPolls <- list()
usaTodayPolls <- list()

length( apPolls ) <- nrow( featureTable2013 )
length( usaTodayPolls ) <- nrow( featureTable2013 )

for( i in 1:length( apPolls ) )
  {
  apPolls[[i]] <- rep( 0, length( weeks ) )
  usaTodayPolls[[i]] <- rep( 0, length( weeks ) )
  }

cat( "Reading USA Today and AP poll data\n" )

for( y in years )
  {
  cat( "  Season ", y, ": week ", sep = "" )
  for( w in weeks )
    {
    cat( w )
    pollUrl <- paste0( "http://espn.go.com/mens-college-basketball/rankings/_/", "year/", y, "/week/", w, "/seasontype/2" )
    pageTables <- readHTMLTable( pollUrl )

    # AP Poll

    pageTables[[1]] <- subset( pageTables[[1]], select = c( "V2", "V4" ) )
    pageTables[[1]] <- pageTables[[1]][3:length( pageTables[[1]][,1]),]
    colnames( pageTables[[1]] ) <- c( "Team", "POINTS" )

    pageTables[[1]] <- data.frame( lapply( pageTables[[1]], as.character ), stringsAsFactors = FALSE )

    for( i in 1:length( pageTables[[1]]$Team ) )
      {
      q <- strsplit( as.character( pageTables[[1]]$Team[i] ), "\\(" )[[1]]
      if( length( q ) > 1 )
        {
        r <- strsplit( q[2], "\\)" )[[1]][1]
        if( !is.na( as.numeric( r ) ) )
          {
          pageTables[[1]]$Team[i] <- gsub( "\\s", "", q[1] )
          }
        }
      idx <- which( featureTable2013$Team == pageTables[[1]]$Team[i] )
      if( length( idx ) > 0 )
        {
        apPolls[[idx]][w] <- as.numeric( gsub( ",", "", pageTables[[1]]$POINTS[i] ) )
        }
      }

    # USA Today poll

    pageTables[[2]] <- subset( pageTables[[2]], select = c( "V2", "V4" ) )
    pageTables[[2]] <- pageTables[[2]][3:length( pageTables[[2]][,1]),]
    colnames( pageTables[[2]] ) <- c( "Team", "POINTS" )

    pageTables[[2]] <- data.frame( lapply( pageTables[[2]], as.character ), stringsAsFactors = FALSE )

    for( i in 1:length( pageTables[[2]]$Team ) )
      {
      q <- strsplit( as.character( pageTables[[2]]$Team[i] ), "\\(" )[[1]]
      if( length( q ) > 1 )
        {
        r <- strsplit( q[2], "\\)" )[[1]][1]
        if( !is.na( as.numeric( r ) ) )
          {
          pageTables[[2]]$Team[i] <- gsub( "\\s", "", q[1] )
          }
        }
      idx <- which( featureTable2013$Team == pageTables[[2]]$Team[i] )
      if( length( idx ) > 0 )
        {
        usaTodayPolls[[idx]][w] <- as.numeric( gsub( ",", "", pageTables[[2]]$POINTS[i] ) )
        }
      }

    if( w == weeks[length(weeks)] )
      {
      cat( "\n", sep = "" )
      }
    else
      {
      cat( ", ", sep = "" )
      }
    }
  }

pollHeaders <- c( "AP.AverageVotes.1to9weeks", "AP.AverageVotes.10to18weeks", "AP.Trend.13to18weeks",
                  "USAToday.AverageVotes.1to9weeks", "USAToday.AverageVotes.10to18weeks", "USAToday.Trend.13to18weeks"
                )
featureTable2013[pollHeaders] <- 0
for( i in 1:nrow( featureTable2013 ) )
  {
  featureTable2013$AP.AverageVotes.1to9weeks[i] <- mean( apPolls[[i]][1:9], na.rm = TRUE )
  featureTable2013$AP.AverageVotes.10to18weeks[i] <- mean( apPolls[[i]][10:18], na.rm = TRUE )
  featureTable2013$USAToday.AverageVotes.1to9weeks[i] <- mean( usaTodayPolls[[i]][1:9], na.rm = TRUE )
  featureTable2013$USAToday.AverageVotes.10to18weeks[i] <- mean( usaTodayPolls[[i]][10:18], na.rm = TRUE )

  spl <- smooth.spline( weeks, apPolls[[i]], spar = 0.5 )
  spl.deriv <- predict( spl, weeks, deriv = 1 )
  featureTable2013$AP.Trend.13to18weeks[i] <- mean( spl.deriv$y[13:18], na.rm = TRUE )

  spl <- smooth.spline( weeks, usaTodayPolls[[i]], spar = 0.5 )
  spl.deriv <- predict( spl, weeks, deriv = 1 )
  featureTable2013$USAToday.Trend.13to18weeks[i] <- mean( spl.deriv$y[13:18], na.rm = TRUE )
  }


##################################################################
#
# Build model
#
##################################################################

##################################################################
#
# Run leave-one-out cross validation
#   1. Use each season as testing data with the remaining entries
#      as training data in a random decision forest
#
##################################################################


numberOfTrees <- 1000

cat( "Constructing random forests for each season using n = ", numberOfTrees, " trees.\n", sep = "" )

cat( "  Season 2013...", sep = "" )
testingData <- subset( featureTable2013, select = -c( Team ) )
trainingData <- subset( featureTable, select = -c( Team, Season ) )

forest <- randomForest( as.formula( "PostSeasonTotalPoints ~ ." ), data = trainingData,
                        na.action = na.omit, ntree = numberOfTrees, importance = FALSE,
                        type = regression )

# forestImp <- importance( forest )
# forestImp.df <- data.frame( Statistic = names( forestImp[,1] ), Importance = as.numeric( forestImp[,1] )  )
# forestImp.df <- forestImp.df[order( forestImp.df$Importance ),]

testResponse <- predict( forest, testingData )

results <- data.frame( Teams = teams2013, PostSeasonPoints = testResponse )
results[order( results$PostSeasonPoints ),]

cat( "  \n", sep = "" )
