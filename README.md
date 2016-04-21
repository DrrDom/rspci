SPCI software calculates fragments contributions to an investigated property by means of interpretation of QSAR models. This package helps to analyze and visualize contributions.

# How to install

First install `devtools` package
```
instal.packages("devtools")
```
Then run from the R console
```
devtools::install_github("DrrDom/rspci")
```
Use the same command to update the package if necessary.

# How to use

#### Load package and prepare data before visualization
```
library(rspci)
```

load text file with contributions (normally this file is located in a dir of modeled property and names default_frag_contributions.txt or auto_frag_contributions.txt, etc)
```
file_name <- system.file("extdata", "BBB_frag_contributions.txt", package = "rspci")
d <- load_data(file_name)
```

add full names containing info about fragments occurencies
```
d <- add_full_names(d)
```

filter fragments by count to discard rarely occured fragments (optional); there are other filters which may be applied
```
d <- filter_by_frags_count(d)
```

reorder data according to specific model and contribution for nice visualization and consistent order of fragments on different plots
```
d <- reorder_data(d, "consensus", "overall")
```

add significance levels to data
```
d <- add_signif(d)
```


#### Visualization of overall contributions
keep only overall contributions for visualization
```
d1 <- filter_by_prop_names(d, "overall")
```

boxplot
```
plot_contrib(d1)
```
barplot
```
plot_contrib(d1, plot_type = "barplot")
```
flipped boxplot with significance levels
```
plot_contrib(d1, flip = FALSE, show_sign_text = "ptext")
```
barpplot with significance levels
```
plot_contrib(d1, plot_type = "barplot", show_sign_text = "ptext")
```


#### Visualization of physico-chemical contributions
keep only physico-chemical contributions
```
d2 <- filter_by_prop_names(d, remove_prop_names = "overall")
```
barplot
```
plot_contrib(d2, plot_type = "barplot")
```
not fliped barplot
```
plot_contrib(d2, plot_type = "barplot", flip = FALSE)
```
not fliped barplot with siginificance levels
```
plot_contrib(d2, plot_type = "barplot", flip = FALSE, show_sign_text = "ptext")
```


#### Create summary table
These commands return summary table of the consensus model sorted by median overall contributions of the fragments
```
library(dplyr)  # to use %>% function

d <- add_full_names(d)  # required, as full_name column will be used for grouping
d <- filter_by_frags_count(d)  # optional, to remove rarely occured fragments

df <- d %>%
  filter(Property == "overall", Model == "consensus") %>%
  group_by(full_name) %>%
  summarise(median = median(Contribution)) %>%
  arrange(desc(median))
```
