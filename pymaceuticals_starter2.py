#!/usr/bin/env python
# coding: utf-8

# ## Observations and Insights 

# 

# In[97]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np


# In[98]:


# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
combined_study_data = pd.merge(mouse_metadata, study_results, on="Mouse ID", how="outer")

# Display the data table for preview
combined_study_data.head()


# In[99]:


# Checking the number of mice.
mice=combined_study_data["Mouse ID"].value_counts()
number_of_mice=len(mice)
number_of_mice


# In[100]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mice = combined_study_data.loc[combined_study_data.duplicated(subset=['Mouse ID', 'Timepoint',]),'Mouse ID'].unique()


# In[101]:


# Optional: Get all the data for the duplicate mouse ID. 
all_duplicate_mouse_id=pd.DataFrame(duplicate_mice)
all_duplicate_mouse_id


# In[102]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_df = combined_study_data[combined_study_data['Mouse ID'].isin(duplicate_mice)==False]


# In[103]:


# Checking the number of mice in the clean DataFrame.
clean_mice=clean_df["Mouse ID"].value_counts()
clean_number_of_mice=len(clean_mice)
clean_number_of_mice


# ## Summary Statistics

# In[104]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# mean, median, variance, standard deviation, and SEM of the tumor volume. 


# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 


# Assemble the resulting series into a single summary dataframe.


# In[105]:


# mean

drug_regimen_mean = clean_df.groupby('Drug Regimen').mean()["Tumor Volume (mm3)"]
drug_regimen_mean


# In[106]:


# median

drug_regimen_median = clean_df.groupby("Drug Regimen").median()["Tumor Volume (mm3)"]
drug_regimen_median


# In[107]:


# variance

drug_regimen_variance = clean_df.groupby("Drug Regimen").var()["Tumor Volume (mm3)"]
drug_regimen_variance


# In[108]:


# standard deviation

drug_regimen_standard_deviation = clean_df.groupby("Drug Regimen").std()["Tumor Volume (mm3)"]
drug_regimen_standard_deviation


# In[109]:


# SEM

drug_regimen_sem = clean_df.groupby("Drug Regimen").sem()["Tumor Volume (mm3)"]
drug_regimen_sem


# In[110]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

summary_statistics_table = pd.DataFrame({"Mean": drug_regimen_mean, "Median": drug_regimen_median, "Variance": drug_regimen_variance, "Standard Deviation": drug_regimen_standard_deviation, "SEM": drug_regimen_sem})

summary_statistics_table


# In[111]:


# Using the aggregation method, produce the same summary statistics in a single line

single_group_by = clean_df.groupby('Drug Regimen')
summary_statistics_table_2 = single_group_by.agg(['mean', 'median', 'variance', 'standard_deviation', 'SEM'])["Tumor Volume (mm3)"]
summary_statistics_table_2


# ## Bar and Pie Charts

# In[112]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.

count_mice_per_treatment = combined_study_data.groupby(["Drug Regimen"]).count()["Mouse ID"]

plot_pandas = count_mice_per_treatment.plot.bar(figsize=(16,10), color='b', fontsize = 16)
count_mice_per_treatment
plt.xlabel("Drug Regimen", fontsize=12)
plt.ylabel("Treatment Cohort Count", fontsize=12)
plt.title("Number of Mice per Treatment Cohort", fontsize=16)

plt.savefig(bbox_inches="tight")
plt.tight_layout()
plt.show()

count_mice_per_treatment


# In[113]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.

drug_regimen = ["Capomulin", "Ceftamin", "Infubinol", "Ketapril", "Naftisol", "Placebo", "Propriva", "Ramicane", "Stelasyn", "Zoniferol"]
members = [225, 175, 175, 185, 180, 178, 165, 225, 185, 190]


# In[117]:


x_axis = np.arrange(0, len(drug_regimen))
# tick_locations = [x for x in x_axis]
# tick_locations = [value for value in x_axis]
tick_locations = []
for x in x_axis:
    tick_locations.append(x)
    
plt.title("Drug Regimen Count")
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Members")

plt.xlim(-0.75, len(drug_regimen)-0.25)
plt.ylim(0, max(members) +5)

plt.bar(x_axis, drug_regimen, facecolor="red", alpha=1, align="center")
plt.xtick(tick_locations, members)
plt.show()


# In[118]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas
groupby_gender = combined_study_data.groupby(["Mouse ID", "Sex"])
groupby_gender
gender_df = pd.DataFrame(groupby_gender.size())

mouse_gender = pd.DataFrame(gender_df.groupby(["Sex"]).count())
mouse_gender.columns - ["Total Count"]

mouse_gender["Percentage of Sex"] = (100*(mouse_gender["Total Count"]/mouse_gender["Total Count"].sum()))

mouse_gender["Percentage of Sex"] = mouse_gender["Percentage of Sex"]

mouse_gender



# In[115]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
labels = ["Female", "Male"]


# ## Quartiles, Outliers and Boxplots

# In[124]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin

capomulin_df = combined_study_data.loc[combined_study_data["Drug Regimen"] == "capomulin",:]

# Start by getting the last (greatest) timepoint for each mouse
capomulin_last = capomulin_df.groupby("Mouse ID").max()["Timepoint"]
capomulin_vol = pd.DataFrame(capomulin_last)
capomulin_merge = pd.merge(capomulin_vol, combined_study_data, on=("Mouse ID","Timepoint"),how="left")
capomulin_merge.head()

# Merge this group df with the original dataframe to get the tumor volume at the last timepoint


# In[125]:


# Ramicane
ramicane_df = combined_study_data.loc[combined_study_data["Drug Regimen"] == "ramicane", :]

ramicane_last = ramicane_df.groupby("Mouse ID").max()["Timepoint"]
ramicane_vol = pd.DataFrame(ramicane_last)
ramicane_merge = pd.merge(ramicane_vol, combined_study_data, on=("Mouse ID","Timepoint"),how="left")
ramicane_merge.head()


# In[126]:


# Infubinol
infubinol_df = combined_study_data.loc[combined_study_data["Drug Regimen"] == "infubinol", :]

infubinol_last = infubinol_df.groupby("Mouse ID").max()["Timepoint"]
infubinol_vol = pd.DataFrame(infubinol_last)
infubinol_merge = pd.merge(infubinol_vol, combined_study_data, on=("Mouse ID","Timepoint"),how="left")
infubinol_merge.head()


# In[127]:


# Ceftamin
ceftamin_df = combined_study_data.loc[combined_study_data["Drug Regimen"] == "ceftamin", :]

ceftamin_last = ceftamin_df.groupby("Mouse ID").max()["Timepoint"]
ceftamin_vol = pd.DataFrame(ceftamin_last)
ceftamin_merge = pd.merge(ceftamin_vol, combined_study_data, on=("Mouse ID","Timepoint"),how="left")
ceftamin_merge.head()


# In[128]:


# Put treatments into a list for for loop (and later for plot labels)


# Create empty list to fill with tumor vol data (for plotting)


# Calculate the IQR and quantitatively determine if there are any potential outliers. 

    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    
    
    # add subset 
    
    
    # Determine outliers using upper and lower bounds
    


# In[15]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest


# ## Line and Scatter Plots

# In[16]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin


# In[17]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen


# ## Correlation and Regression

# In[18]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen


# In[ ]:




