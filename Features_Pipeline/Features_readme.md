# MATLAB Feature Extraction and Analysis 

## How to run using example data
1. Download the example from here: https://drive.google.com/file/d/1ifZrDXY1Ib7znB6O0f8xXPZ0iNo1Hd72/view?usp=drive_link
2. Extract it in this current folder
3. In MATLAB, add to path the _data_, _Nature_, and _spider_plot-master_ directories
4. In Correlations Runners folder, open DefaultRunner.m
5. Run it cell by cell. The results are saved by default to _C:\CellInsights\Results\Example_
   6. This results path can be modified at the start of DefaultRunner.m

## How to create new dataset from the X-Y positions (i.e. from the Tracker output)
1. Important - This below is currently manually done. Automated process is to be developed/ 
2. Import the output.csv X-Y trajectories file from its location into MATLAB as _Column Vectors_
   2. Default location is: C:\CellInsights\Tracking
3. Open Main_Part1.m 
   4. Change the default path to your path at the start of this file
   5. Run all cells
   6. For the MEAN SQUARED DISPLACEMENT cell, choose the newly created _Data_CellID_Frame_X_Y.xlsx_ file. 
      7. The dialog also asks for a trajectory time step size. 
         8. For example, in our default research we took images of the cells every 10 minutes, so we use 10 at this point
      9. Save the created _MSD.xlsx_ file (you will be prompt to save), preferably at the default directory
   10. For the  FITTING OF CELL TRAJECTORIES TO THE APRW MODEL, also choose the newly created _Data_CellID_Frame_X_Y.xlsx_ file. 
       11. The dialog also asks for a trajectory time step size (10 is our default).
       12. Save the created _APRW model fit.xlsx_ file (you will be prompt to save), preferably at the default directory
   13. Dont forget to run the last matlab cell
   14. DONT DELETE THE VARIABLES, WE WILL NEED THEM SOON
15. Open and save a new copy of _APRW Model Features - Empty.xlsx_, it can be found in this current foler
16. Fill it with all needed information:
    17. Pp and Sp from the two first columns of sheet _APRW fitting-p_ of _APRW model fit.xlsx_
    18. Pnp and Snp from the two first columns of sheet _APRW fitting-np_ of _APRW model fit.xlsx_
    19. MSD information:
        20. Open MSD.xlsx
        21. copy all data from _individual cell msd sheet_
        22. Create a new sheet and paste the data transposed
        23. Create a new column and Average all cells by columns
        24. Copy this new created column to the MSD part of the _APRW Model Features - Empty.xlsx_
    25. In MATLAB, copy the _CELLTYPE_unique_ variable as a column. Paste it at the _CELLTYPE_ part of the _APRW Model Features - Empty.xlsx_
    26. in MATLAB, copy the _CellID_unique_ variable as a column. Paste it at the _CellId_ part of the _APRW Model Features - Empty.xlsx_
    27. In MATLAB, copy the _CELLID_unique_ variable as a column. Paste it at the _Patient/CELLID2_ part of the _APRW Model Features - Empty.xlsx_
    28. Manually enter the CAP score per patient. For the sake of the example all CAP scores can be 0.
    29. Save the _APRW Model Features - Empty.xlsx_ file
30. Delete all variables in MATLAB.
31. Import the _APRW Model Features - Empty.xlsx_ file into MATLAB as _Column Vectors_
32. Open Main_Part2.m
    33. Change the default path to your path at the start of this file
    34. Run all cells
    35. For the TRJs APRW cell, choose the _APRW model fit.xlsx_ file
        36. The dialog also asks for a trajectory time step size (10 is our default). 
    37. For the TRJs Full cell, choose the _Data_CellID_Frame_X_Y.xlsx_ file
    38. Dont forget to run the last cell
39. In MATLAB, add to path the _Nature_, and _spider_plot-master_ directories
40. In Correlations Runners folder, open DefaultRunner.m
41. Activate the section called "Activate the below if using own data"
    42. Disable the section called "Activate the below if using example data"
43. Run the notebook