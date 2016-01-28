SPCI software calculates fragments contributions to an investigated property by means of interpretation of QSAR models. This package helps to analyze and visualize contributions.

# Help

#### Load package prepare data before visualization
```
library(rspci)
```

load text file with contributions
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
boxplot with siginificance levels
```
plot_contrib(d2, show_sign_text = "ptext")
```
