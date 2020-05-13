## Multi-population Modified L-SHADE for Single Objective Bound Constrained Optimization(mpmL-SHADE)

### About the authors:
 The program mpmL-SHADE is implemented by Yann-Chern Jou, Shuo-Ying Wang, Jia-Fong Yeh, and Tsung-Che Chiang.  
 
 If you have any question about the algorithm or program, please contact us.

### About the algorithm:
 We participated the competitions of Real-World Single Objective Constrained Optimization in CEC Special Session, IEEE WCCI 2020. The paper of our algorithm is published in CEC2020. We extend a previous algorithm mL-SHADE by running the evolutionary process through multiple populations and adding dynamic control of mutation intensity and hyper-parameters. The whole population is partitioned into sub-populations by a random clustering method. Mutation intensity and hyper-parameters are adjusted based on the consumption of fitness function evaluations. Performance of the proposed algorithm is verified by ten benchmark functions in the [CEC2020 Competition on Single Objective Bound Constrained Optimization](https://github.com/P-N-Suganthan/2020-Bound-Constrained-Opt-Benchmark "link").

### About the program:
 We developed our program by C++ (std C++17). The development environment is Visual Studio 2019 under Win10(SDK: 10.0.18363)  
 
 The detail of executation is shown down below:  
   
 * __Step 1.__ Set the program parameter in "problem_list_DXX.ini"  
         Currently setting is CEC2020 benchmark test problem for 4 different dimensions(D = 5, 10, 15, 20).  
         
 * __Step 2.__ Set the experiment list in "experiment_list.ini"  
         
 
 * __Step 3.__ In the 'experiments' folder, we have parameters setting file of mpmL-SHADE for each problem, you can modify the value if you want.
 
 * __Step 4.__ Then, You can execute the 'CEC2020-Bound-Constrained-Opt_mpmL-SHADE.exe' to run the experiments.  
   
   Note: If you also use our log system, please make sure the storage capacity of disk is more than 100GB.  
 
### About the results:
  Please check the "results" folder. The more detail of statistical analysis is in the section IV. of our CEC paper.  
  All of computing results were proccessed by our helping program in the "Analyzing Tools" folders.
