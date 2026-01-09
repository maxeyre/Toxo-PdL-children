# To plot two series vertically aligned with only one labelled x-axis,
# we can remove the axes from the top plot and then plot the two graphs
# together using either the egg package or the cowplot package

library(ggplot2) # 3.2.1
library(dplyr)
library(lubridate)
library(cowplot) # 1.0.0
library(egg) # 0.4.5

#' Create some data to play with. Two time series with the same timestamp.
df <- data.frame(DateTime = ymd("2010-07-01") + c(0:8760) * hours(2), 
                 series1 = rnorm(8761), 
                 series2 = rnorm(8761, 100))

#' Create the two plots.
plot1 <- df %>%
  select(DateTime, series1) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series1), size = 0.5, alpha = 0.75) +
  ylab("Red dots / m") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

plot2 <- df %>%
  select(DateTime, series2) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series2), size = 0.5, alpha = 0.75) +
  ylab("Blue drops / L") +
  theme_minimal() +
  theme(axis.title.x = element_blank())

# Draw the two plot aligned vertically, with the top plot 1/3 of the height
# of the bottom plot
cowplot::plot_grid(plot1, plot2, align = "v", ncol = 1, rel_heights = c(0.25, 0.75))
egg::ggarrange(plot1, plot2, heights = c(0.25, 0.75))



#' When plotting multiple data series that share a common x axis but different y axes,
#' we can just plot each graph separately. This suffers from the drawback that the shared axis will typically
#' not align across graphs due to different plot margins.
#' One easy solution is to reshape2::melt() the data and use ggplot2's facet_grid() mapping. However, there is
#' no way to label individual y axes.
#' facet_grid() and facet_wrap() were designed to plot small multiples, where both x- and y-axis ranges are
#' shared acros all plots in the facetting. While the facet_ calls allow us to use different scales with
#' the \code{scales = "free"} argument, they should not be used this way.
#' A more robust approach is to the grid package grid.draw(), rbind() and ggplotGrob() to create a grid of 
#' individual plots where the plot axes are properly aligned within the grid.
#' Thanks to https://rpubs.com/MarkusLoew/13295 for the grid.arrange() idea.

library(ggplot2)
library(grid)
library(dplyr)
library(lubridate)

#' Create some data to play with. Two time series with the same timestamp.
df <- data.frame(DateTime = ymd("2010-07-01") + c(0:8760) * hours(2), series1 = rnorm(8761), series2 = rnorm(8761, 100))

#' Create the two plots.
plot1 <- df %>%
  select(DateTime, series1) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series1), size = 0.5, alpha = 0.75) +
  ylab("Red dots / m") +
  theme_minimal() +
  theme(axis.title.x = element_blank())

plot2 <- df %>%
  select(DateTime, series2) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series2), size = 0.5, alpha = 0.75) +
  ylab("Blue drops / L") +
  theme_minimal() +
  theme(axis.title.x = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
