library(Lahman)
library(tidyverse)
library(dslabs)
ds_theme_set()
# Teams %>% filter(yearID %in% 1961:2001) %>%
#   mutate(HR_per_game = HR / G, R_per_game = R / G) %>%
#   ggplot(aes(HR_per_game, R_per_game)) +
#   geom_point(alpha = 0.5)
# Teams %>% filter(yearID %in% 1961:2001) %>%
#   mutate(SB_per_game = SB / G, R_per_game = R / G) %>%
#   ggplot(aes(SB_per_game, R_per_game)) +
#   geom_point(alpha = 0.5)
# Teams %>% filter(yearID %in% 1961:2001) %>%
#   mutate(BB_per_game = BB / G, R_per_game = R / G) %>%
#   ggplot(aes(BB_per_game, R_per_game)) +
#   geom_point(alpha = 0.5)
# Teams %>% filter(yearID %in% 1961:2001) %>%
#   mutate(AA_per_game = E / G, R_per_game = R / G) %>%
#   ggplot(aes(AA_per_game, R_per_game)) +
#   geom_point(alpha = 0.5)
# Teams %>% filter(yearID %in% 1961:2001) %>%
#   ggplot(aes(X3B, X2B)) +
#   geom_point(alpha = 0.5)
# my_team <- Teams %>% filter(yearID %in% 1961:2001) %>%
#   mutate(R_per_game = R / G,AB_per_game = AB / G) %>%
#   ggplot(aes(R_per_game,AB_per_game)) +
#   geom_point(alpha = 0.5)
# my_team <- Teams %>% filter(yearID %in% 1961:2001) %>%
#   mutate(R_per_game = R / G, AB_per_game = AB / G)
# cor(my_team$R_per_game,my_team$AB_per_game)
# my_team <- Teams %>% filter(yearID %in% 1961:2001) %>%
#   mutate(W_per_game = W / G, ER_per_game = E / G)
# cor(my_team$W_per_game,my_team$ER_per_game)
# my_team <- Teams %>% filter(yearID %in% 1961:2001) %>%
#   mutate(X2B_per_game = X2B / G, X3B_per_game = X3B / G)
# cor(my_team$X2B_per_game,my_team$X3B_per_game)
# set.seed(1989) #if you are using R 3.5 or earlier
# set.seed(1989, sample.kind="Rounding") #if you are using R 3.6 or later

# 
# female_heights <- GaltonFamilies%>%     
#   filter(gender == "female") %>%     
#   group_by(family) %>%     
#   sample_n(1) %>%     
#   ungroup() %>%     
#   select(mother, childHeight) %>%     
#   rename(daughter = childHeight)
# mu_x <- mean(female_heights$mother)
# s_x <- sd(female_heights$mother)
# mu_y <- mean(female_heights$daughter)
# s_y <- sd(female_heights$daughter)
# r <- cor(female_heights$mother,female_heights$daughter)
# m_1 <-  r * s_y / s_x
# b_1 <- mu_y - m_1*mu_x
# m_1
# b_1
# x <- r*s_x/s_y
# y <- r^2*100
# y
# z <- m_1*60+b_1
# z
# library(HistData)
# data("GaltonFamilies")
# set.seed(1989, sample.kind="Rounding") #if you are using R 3.6 or later
# library(HistData)
# data("GaltonFamilies")
# options(digits = 3)    # report 3 significant digits
# 
# female_heights <- GaltonFamilies %>%     
#   filter(gender == "female") %>%     
#   group_by(family) %>%     
#   sample_n(1) %>%     
#   ungroup() %>%     
#   select(mother, childHeight) %>%     
#   rename(daughter = childHeight)
# fit <- lm(female_heights$mother ~ female_heights$daughter , data = female_heights)
# summary(fit)
# X1 <- fit$coefficients[1]+fit$coefficients[2]*female_heights$daughter[1]
# X1
# female_heights$mother[1]
# library(Lahman)
# bat_02 <- Batting %>% filter(yearID == 2002) %>%
#   mutate(pa = AB + BB, singles = (H - X2B - X3B - HR)/pa, bb = BB/pa) %>%
#   filter(pa >= 100) %>%
#   select(playerID, singles, bb)
# bat_99_01 <- Batting %>% filter(yearID %in% 1999:2001) %>%
#   mutate(pa = AB + BB, singles = (H - X2B - X3B - HR)/pa, bb = BB/pa) %>%
#   filter(pa >= 100) %>%
#   group_by(playerID) %>%
#   summarize(mean_singles = mean(singles), mean_bb = mean(bb))
# sum(bat_99_01$mean_singles > 0.2)
# sum(bat_99_01$mean_bb > 0.2)
# dat <- inner_join(bat_02, bat_99_01)
# cor(dat$singles, dat$mean_singles)
# library(tidyverse)
# library(HistData)
# data("GaltonFamilies")
# set.seed(1) # if you are using R 3.5 or earlier
# set.seed(1, sample.kind = "Rounding") # if you are using R 3.6 or later
# galton <- GaltonFamilies %>%
#   group_by(family, gender) %>%
#   sample_n(1) %>%
#   ungroup() %>% 
#   gather(parent, parentHeight, father:mother) %>%
#   mutate(child = ifelse(gender == "female", "daughter", "son")) %>%
#   unite(pair, c("parent", "child"))
# 
# X <- galton %>%
#   group_by(pair) %>%
#   summarize(cor = cor(parentHeight, childHeight))%>%
#   filter(cor == max(cor))
# library(broom)
# YY <- galton %>%
#   group_by(pair) %>%
#   do(tidy(lm(childHeight ~ parentHeight, data = .), conf.int = TRUE))
# 
# library(broom)
# galton %>%
#   group_by(pair) %>%
#   do(tidy(lm(childHeight ~ parentHeight, data = .), conf.int = TRUE)) %>%
#   filter(term == "parentHeight", pair == "father_daughter") %>%
#   pull(estimate)
# galton %>%
#   group_by(pair) %>%
#   do(tidy(lm(childHeight ~ parentHeight, data = .), conf.int = TRUE)) %>%
#   filter(term == "parentHeight", pair == "mother_son") %>%
#   pull(estimate)
# 
# galton %>%
#   group_by(pair) %>%
#   do(tidy(lm(childHeight ~ parentHeight, data = .), conf.int = TRUE)) 
# library(Lahman)
# library(broom)
# Teams %>%
#   filter(yearID == 1971) %>%
#   lm(R ~ BB + HR, data = .) %>%
#   tidy() %>%
#   filter(term == "HR") %>%
#   pull(estimate)
# 
# Teams %>%
#   filter(yearID == 2018) %>%
#   lm(R ~ BB + HR, data = .) %>%
#   tidy() %>%
#   filter(term == "BB") %>%
#   pull(estimate)
# res <- Teams %>%
#   filter(yearID %in% 1961:2018) %>%
#   group_by(yearID) %>%
#   do(tidy(lm(R ~ BB + HR, data = .))) %>%
#   ungroup() 
# res %>%
#   filter(term == "BB") %>%
#   ggplot(aes(yearID, estimate)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# res <- Teams %>%
#   filter(yearID %in% 1961:2018) %>%
#   group_by(yearID) %>%
#   do(tidy(lm(R ~ BB + HR, data = .))) %>%
#   ungroup() 
# X <- res %>% filter(term == "BB") %>%
#   lm(estimate ~ yearID, data = .) %>%
#   tidy()%>%
# res %>%
#   filter(term == "BB") %>%
#   lm(estimate ~ yearID, data = .) %>%
#   tidy() %>%
#   filter(term == "yearID") %>%
#   pull(p.value)
# res %>%
#   filter(term == "BB") %>%
#   lm(estimate ~ yearID, data = .) %>%
#   tidy() %>%
#   filter(term == "yearID") %>%
# #   pull(estimate)
# library(tidyverse)
# library(broom)
# library(Lahman)
# Teams_small <- Teams %>% 
#   filter(yearID %in% 1961:2001) %>% 
#   mutate(avg_attendance = attendance/G , RG = R/G , HRG = HR/G)
# head(Teams_small)
# Teams_small %>% 
#   lm(avg_attendance ~ RG , data = .) %>% 
#   tidy()%>%
#   pull(estimate)
# Teams_small %>% 
#   lm(avg_attendance ~ HRG , data = .) %>% 
#   tidy()%>%
#   pull(estimate)
# Teams_small %>% 
#   lm(avg_attendance ~ W , data = .) %>% 
#   tidy()%>%
#   pull(estimate)
# Teams_small %>% 
#   lm(avg_attendance ~ yearID , data = .) %>% 
#   tidy()%>%
#   pull(estimate)
# Teams_small %>% summarise(cor= cor(W , HRG )) %>%
#   pull(cor)
# dat <- Teams_small %>%
#   mutate(W_strata = round(W/10)) %>%
#   filter(W_strata >= 5 & W_strata <= 10) %>%
#   group_by(W_strata) %>%
#   do(tidy(lm(avg_attendance ~ RG , data=.), int.cof = TRUE)) %>%
#   filter(term == "RG")
# dat
# Teams_small %>% 
#   lm(avg_attendance ~ RG +HRG + W + yearID , data = .) %>% 
#   tidy()%>%
#   pull(estimate)
# fit <- Teams_small %>% 
#   mutate(R_per_game = R/G,
#          HR_per_game = HR/G) %>%
#   lm(avg_attendance ~ R_per_game + HR_per_game + W + yearID, data = .)
# tidy(fit) %>%
#   filter(term == "R_per_game") %>%
#   pull(estimate)
# dat1 <- Teams_small %>%
#   mutate(W_strata = round(W/10)) %>%
#   filter(W_strata >= 5 & W_strata <= 10) %>%
#   group_by(W_strata) %>%
#   do(tidy(lm(avg_attendance ~ HRG , data=.), int.cof = TRUE)) %>%
#   filter(term == "HRG")
# 
# predict(fit, data.frame(R_per_game = 5, HR_per_game = 1.2, W = 80, yearID = 2002))
# 
# fit <- Teams%>% 
#   mutate(R_per_game = R/G,
#          HR_per_game = HR/G) %>%
#   lm(avg_attendance ~ R_per_game + HR_per_game + W + yearID, data = .)
# tidy(fit) 
# x <- predict(fit, data.frame(yearID = 2002))
# y <- Teams
# newdata <- Teams %>%
#   filter(yearID == 2002) %>%
#   mutate(avg_attendance = attendance/G,
#          R_per_game = R/G,
#          HR_per_game = HR/G)
# preds <- predict(fit, newdata)
# cor(preds, newdata$avg_attendance)
# N <- 25
# g <- 1000000
# sim_data <- tibble(group = rep(1:g, each = N), x = rnorm(N * g), y = rnorm(N * g))
# 
# # calculate correlation between X,Y for each group
# res <- sim_data %>% 
#   group_by(group) %>% 
#   summarize(r = cor(x, y)) %>% 
#   arrange(desc(r))
# res
# 
# # plot points from the group with maximum correlation
# sim_data %>% filter(group == res$group[which.max(res$r)]) %>%
#   ggplot(aes(x, y)) +
#   geom_point() + 
#   geom_smooth(method = "lm")
# 
# # histogram of correlation in Monte Carlo simulations
# res %>% ggplot(aes(x=r)) + geom_histogram(binwidth = 0.1, color = "black")
# 
# # linear regression on group with maximum correlation
# library(broom)
# sim_data %>% 
#   filter(group == res$group[which.max(res$r)]) %>%
#   do(tidy(lm(y ~ x, data = .)))
# 
# # simulate independent X, Y and standardize all except entry 23
# set.seed(1985)
# x <- rnorm(100,100,1)
# y <- rnorm(100,84,1)
# x[-23] <- scale(x[-23])
# y[-23] <- scale(y[-23])
# 
# # plot shows the outlier
# qplot(x, y, alpha = 0.5)
# 
# # outlier makes it appear there is correlation
# cor(x,y)
# cor(x[-23], y[-23])
# 
# # use rank instead
# qplot(rank(x), rank(y))
# cor(rank(x), rank(y))
# 
# # Spearman correlation with cor function
# cor(x, y, method = "spearman")
library(dslabs)
data("research_funding_rates")
research_funding_rates
two_by_two <- research_funding_rates %>% 
  select(-discipline) %>% 
  summarize_all(funs(sum)) %>%
  summarize(yes_men = awards_men, 
            no_men = applications_men - awards_men, 
            yes_women = awards_women, 
            no_women = applications_women - awards_women) %>%
  gather %>%
  separate(key, c("awarded", "gender")) %>%
  spread(gender, value)
two_by_two %>% 
  mutate(men = round(men/sum(men)*100, 1), women = round(women/sum(women)*100, 1)) %>%
  filter(awarded == "yes") %>%
  pull(men)
dat <- research_funding_rates %>% 
  mutate(discipline = reorder(discipline, success_rates_total)) %>%
  rename(success_total = success_rates_total,
         success_men = success_rates_men,
         success_women = success_rates_women) %>%
  gather(key, value, -discipline) %>%
  separate(key, c("type", "gender")) %>%
  spread(type, value) %>%
  filter(gender != "total")
dat %>% 
  ggplot(aes(discipline, success, size = applications, color = gender)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_point()

dat %>% 
  ggplot(aes(discipline, success, size = applications, color = gender)) + 
    geom_point()