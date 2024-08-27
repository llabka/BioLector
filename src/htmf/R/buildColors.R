buildColors = function(
    wellTypes,
    wellColors = NULL, # rainbow(length(levels(as.factor(wellTypes))))
    colorFcn = rainbow
)
{
   n <- length(wellTypes)
   colorFcts <- as.factor(wellTypes)
   k <- length(levels(colorFcts))

   if( is.null(wellColors) ) {
      wellColors <- c("black", "red")
      if(k > 2)
         wellColors <- colorFcn(k)
   }

   levels(colorFcts) <- wellColors
   colorList <- as.character(colorFcts)
   names(colorList) <- wellTypes
   return(colorList)
}
