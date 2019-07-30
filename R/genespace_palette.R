#' @title Palettes for GENESPACE
#'
#' @description
#' \code{genespace_pal} A set of color palettes that are colorblind safe
#' and work well with genespace plotting.
#'
#' @param n numeric, length 1, the number of colors to return
#' @param name charcter, length 1, name of palette, options are 'all', 'div' (diverging),
#' 'dark', 'medium' and 'light'.
#' @param choose Numeric, length equal to n. Specify the vector of colors from
#' full palette to choose.
#' @param maroons character, length three, the maroon colors in the palette.
#' @param reds character, length three, the maroon colors in the palette.
#' @param oranges character, length three, the maroon colors in the palette.
#' @param yellows character, length three, the maroon colors in the palette.
#' @param ygreens character, length three, the maroon colors in the palette.
#' @param greens character, length three, the maroon colors in the palette.
#' @param cyans character, length three, the maroon colors in the palette.
#' @param lblues character, length three, the maroon colors in the palette.
#' @param blues character, length three, the maroon colors in the palette.
#' @param purples character, length three, the maroon colors in the palette.
#' @param magentas character, length three, the maroon colors in the palette.
#' @param alpha numeric (0..1) indicating the level of opacity
#' @param plotit logical, should a plot be made?
#' @param ... Additional arguments passed to cull_syntenicBlast
#' @details None yet

#' @return A vector of colors of length n
#'
#' @examples
#' \dontrun{
#' genespace_pal(n = 16, name = "all")
#' }
#' @importFrom grDevices colorRampPalette
#' @export
genespace_pal <- function(n,
                          name,
                          choose = NULL,
                          maroons = c("#640000", "#A03F43","#C88C8C"),
                          reds = c("#b71600", "red", "#FFAAAA"),
                          oranges = c("#f95c00", "#ff7f15", "#FFB469"),
                          yellows = c("#cd950c", "#ffc800", "#ffee7a"),
                          ygreens = c("#68a200", "#95de2b", "#cef375"),
                          greens = c("#207910", "#00C00C", "#92f690"),
                          cyans = c("#007266", "#00d0b9", "#99ffe1"),
                          lblues = c("#00446E", "#10B2FF", "#9AE4FF"),
                          blues = c("#00098C", "#3947ff", "#68a4ff"),
                          purples = c("#5200B6", "#9258FF", "#C3B0FF"),
                          magentas = c("#960096", "#ec40ec", "#ff98ff"),
                          alpha = 1,
                          return.function = F,
                          plotit = T){

  coll <- list(maroons, reds, oranges,
               yellows, ygreens, greens,
               cyans, lblues, blues,
               purples, magentas)

  darks = sapply(coll, function(x) x[1])
  meds = sapply(coll, function(x) x[2])
  lights = sapply(coll, function(x) x[3])

  if(name == "dark"){
    pal <- c(darks[1:4],
             meds[7],
             darks[7:8],
             meds[9:10],
             darks[11])
  }
  if(name == "medium"){
    pal <- c(darks[1],
             meds[1:3],
             meds[4],
             darks[6],
             meds[7:11])
  }
  if(name == "light"){
    pal <- c(lights[1:2],
             meds[3:4],
             lights[5:6],
             meds[7],
             lights[8:9],
             lights[11])
  }
  if(name == "div"){
    pal <- with(coll, c(
      blues[2], maroons[1], oranges[2],
      lblues[3], reds[1], yellows[2],
      purples[3], blues[1], cyans[2],
      magentas[2], maroons[3], blues[3],
      oranges[3], greens[3], purples[2]))
  }
  if(name == "all"){
    pal <- unlist(coll)
  }

  if(n > length(pal)){
    warning("maximum number of colors for palette ", name, " is ",length(pal),"\n",
            "\t... interpolating colors. Colorblind safe not guarunteed.\n")
    pal.fun <- colorRampPalette(pal)
    pal <- pal.fun(n)
  }
  if(n < length(pal)){
    if(is.null(choose)){
      choose <- seq(from = 1,
                    to = length(pal),
                    length.out = n)
    }
    pal <- pal[choose]
  }
  pal <- add_alpha(pal,
                   alpha = alpha)
  if(plotit){
    co <- unlist(coll)
    if(length(co) > n){
      x1 <- 1:length(co)
      x2 <- seq(1,
                length(co),
                length.out = n)
    }else{
      x1 <- seq(from = 1,
                to = n,
                length.out = length(co))
      x2 <- 1:n
    }
    xs <- c(x1, x2)
    ys <- rep(c(1,2),
              c(length(x1), length(x2)))
    cols <- c(co, pal)
    cx <- ifelse(n > 50, 2,
                 ifelse(n > 20, 4, 6))
    plot(xs, ys,
         col = cols,
         cex = cx,
         pch = 15,
         main = paste("full (bottom) vs.", name,"(top) palettes"),
         bty = "n",
         axes = F,
         xlab = "color order",
         ylab = "")
  }
  if(!return.function){
    return(pal)
  }else{
    return(colorRampPalette(pal))
  }

}


