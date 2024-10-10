\ Model Maximize_Kidney_Donations
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  cycle[9,47,31,9] + cycle[31,40,37,31] + cycle[37,31,4,37]
   + cycle[41,37,31,41] + cycle[37,31,41,37] + cycle[40,37,31,40]
   + cycle[9,34,24,9] + cycle[31,25,37,31] + cycle[31,33,37,31]
   + cycle[24,9,34,24] + cycle[5,37,31,5] + cycle[37,31,25,37]
   + cycle[37,31,5,37] + cycle[47,31,9,47] + cycle[31,5,37,31]
   + cycle[34,24,9,34] + cycle[37,31,40,37] + cycle[31,41,37,31]
   + cycle[33,37,31,33] + cycle[25,37,31,25] + cycle[31,9,47,31]
   + cycle[31,4,37,31] + cycle[37,31,33,37] + cycle[4,37,31,4]
Subject To
 Node_Constraint_33: cycle[31,33,37,31] + cycle[33,37,31,33]
   + cycle[37,31,33,37] <= 1
 Node_Constraint_34: cycle[9,34,24,9] + cycle[24,9,34,24]
   + cycle[34,24,9,34] <= 1
 Node_Constraint_4: cycle[37,31,4,37] + cycle[31,4,37,31]
   + cycle[4,37,31,4] <= 1
 Node_Constraint_37: cycle[31,40,37,31] + cycle[37,31,4,37]
   + cycle[41,37,31,41] + cycle[37,31,41,37] + cycle[40,37,31,40]
   + cycle[31,25,37,31] + cycle[31,33,37,31] + cycle[5,37,31,5]
   + cycle[37,31,25,37] + cycle[37,31,5,37] + cycle[31,5,37,31]
   + cycle[37,31,40,37] + cycle[31,41,37,31] + cycle[33,37,31,33]
   + cycle[25,37,31,25] + cycle[31,4,37,31] + cycle[37,31,33,37]
   + cycle[4,37,31,4] <= 1
 Node_Constraint_5: cycle[5,37,31,5] + cycle[37,31,5,37]
   + cycle[31,5,37,31] <= 1
 Node_Constraint_40: cycle[31,40,37,31] + cycle[40,37,31,40]
   + cycle[37,31,40,37] <= 1
 Node_Constraint_9: cycle[9,47,31,9] + cycle[9,34,24,9] + cycle[24,9,34,24]
   + cycle[47,31,9,47] + cycle[34,24,9,34] + cycle[31,9,47,31] <= 1
 Node_Constraint_41: cycle[41,37,31,41] + cycle[37,31,41,37]
   + cycle[31,41,37,31] <= 1
 Node_Constraint_47: cycle[9,47,31,9] + cycle[47,31,9,47]
   + cycle[31,9,47,31] <= 1
 Node_Constraint_24: cycle[9,34,24,9] + cycle[24,9,34,24]
   + cycle[34,24,9,34] <= 1
 Node_Constraint_25: cycle[31,25,37,31] + cycle[37,31,25,37]
   + cycle[25,37,31,25] <= 1
 Node_Constraint_31: cycle[9,47,31,9] + cycle[31,40,37,31]
   + cycle[37,31,4,37] + cycle[41,37,31,41] + cycle[37,31,41,37]
   + cycle[40,37,31,40] + cycle[31,25,37,31] + cycle[31,33,37,31]
   + cycle[5,37,31,5] + cycle[37,31,25,37] + cycle[37,31,5,37]
   + cycle[47,31,9,47] + cycle[31,5,37,31] + cycle[37,31,40,37]
   + cycle[31,41,37,31] + cycle[33,37,31,33] + cycle[25,37,31,25]
   + cycle[31,9,47,31] + cycle[31,4,37,31] + cycle[37,31,33,37]
   + cycle[4,37,31,4] <= 1
Bounds
Binaries
 cycle[9,47,31,9] cycle[31,40,37,31] cycle[37,31,4,37] cycle[41,37,31,41]
 cycle[37,31,41,37] cycle[40,37,31,40] cycle[9,34,24,9] cycle[31,25,37,31]
 cycle[31,33,37,31] cycle[24,9,34,24] cycle[5,37,31,5] cycle[37,31,25,37]
 cycle[37,31,5,37] cycle[47,31,9,47] cycle[31,5,37,31] cycle[34,24,9,34]
 cycle[37,31,40,37] cycle[31,41,37,31] cycle[33,37,31,33]
 cycle[25,37,31,25] cycle[31,9,47,31] cycle[31,4,37,31] cycle[37,31,33,37]
 cycle[4,37,31,4]
End
