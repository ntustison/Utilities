library( ggplot2 )
library( party )

# create two well-separate gaussians
mystd <- 1.0

x1 <- rnorm( 100, mean = 1, sd = mystd )
y1 <- rnorm( 100, mean = 1, sd = mystd )

x2 <- rnorm( 100, mean = 3, sd = mystd )
y2 <- rnorm( 100, mean = 4, sd = mystd )

l <- t( cbind( t( rep( 1, 100 ) ), t( rep( 2, 100 ) ) ) );
myData <- data.frame( cbind( t( cbind( t( x1 ), t( x2 ) ) ), t( cbind( t( y1 ), t( y2 ) ) ), as.factor( l ) ) );
colnames( myData ) <- c( "x", "y", "label" );

myData <- transform( myData, label = factor( label ) );

myPlot <- ggplot( myData, aes( x = x, y = y, group = label ) ) +
          geom_point( data = myData, aes( colour = label, shape = label ), size = 3 )
ggsave( filename = "test.pdf", plot = myPlot, width = 8, height = 6, units = 'in' )

ind <- sample.int( 2, nrow( myData ), replace = TRUE, prob = c( 0.7, 0.3 ) )

myData.train <- myData[ind == 1,]
myData.evaluate <- myData[ind == 2,]


model <- cforest( label ~ x + y, data = myData.train,
  controls = cforest_unbiased( ntree = 100, mtry = 1 ) )

myData.evaluate$prediction <- predict( model, newdata = myData.evaluate )
myData.evaluate$correct <- myData.evaluate$prediction == myData.evaluate$label

response <- treeresponse( model, myData.evaluate )

probs <- matrix( unlist( response, use.names = FALSE ), ncol = 2, byrow = TRUE )

for( j in 1:nrow( myData.evaluate ) )
  {
  myData.evaluate$probabilities[j] <- unlist( response[j] )[as.numeric( myData.evaluate$label[j] )]
  }

varimp( model )
# > test.eva <- as.data.frame( t( c( 2, 2 ) ) )
# > colnames( test.eva )<- c( "x", "y")
# > test.res <- treeresponse( model, test.eva )

#
# myData.evaluate$probabilities1 <- 1 - unlist( treeresponse( model, myData.evaluate ), use.names = F )[seq(1,nrow( myData.evaluate ) * 2, 2)]

myPlot <- ggplot( myData, aes( x = x, y = y, group = label ) ) +
          geom_point( data = myData.train, aes( colour = label, shape = label ), size = 3 ) +
          geom_point( data = myData.evaluate, aes( size = myData.evaluate$probabilities, shape = label ), colour = 'black' )

ggsave( filename = "test.pdf", plot = myPlot, width = 8, height = 6, units = 'in' )



#   thickPlot <- ggplot( plotData, aes( x = Age, y = Thickness, group = Gender ) ) +
#                geom_point( data = plotData, aes( colour = Gender, shape = Gender ), size = 3 ) +
#                scale_x_continuous( "Age (years)", breaks = seq( 20, 90, by = 10 ), labels = seq( 20, 90, by = 10 ), limits = c( 20, 90 ) ) +
#                scale_y_continuous( "Thickness (mm)", breaks = seq( 0, 5, by = 1 ), labels = seq( 0, 5, by = 1 ), limits = c( 0, 5 ) ) +
#                scale_colour_manual( values = c( "navyblue", "darkred" ), breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
#                scale_shape_manual( values = c( 18, 16 ), breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
# #                scale_shape_discrete( breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
#                theme( legend.justification = c( 0, 0 ), legend.position = c( 0, 0 ) ) +
#                ggtitle( paste( "Cortical thickness (", cortical_labels[i], ")", sep = "" ) )
#   ggsave( filename = paste( "yylabel", i, "_results.pdf", sep = "" ), plot = thickPlot, width = 8, height = 6, units = 'in' )
