# ViralKnots

1. First, download the fasta file containing the RNA sequence of the viral genome you wish to run ViralKnots on. The file should be in fasta format and the name of the file should end with .fasta
    - Here is an example of the proper formatting of a fasta file; the file should contain two lines, the first with a greater-than symbol (>) followed by the name of the virus, and the second line containing the viral sequence
    
    ```
    >name_of_virus
    AUGUCUGUCUAUAUCUGUA....
    ```

2. ssh onto sherlock by typing the following command into your terminal: 
    
    `ssh <sunetid>@login.sherlock.stanford.edu -Y`

3. If you do not already have conda installed on sherlock, install it using the following steps:
    - navigate to your home folder and install anaconda there using the following commands (make sure to type 'yes' and accept anything the conda installer requests):

    ```
    cd $HOME:./
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh
    ```

    - now close and reopen your terminal to allow the changes to take effect
    - install biopython & other necessary packages

    ```
    conda install biopython
    conda install scipy
    conda install pandas
    ```

4. in order to use arnie, you must also add the following lines to your bashrc file
    - first, open your bashrc

    `nano ~/.bashrc`
    
    - now, add the following lines to the bottom of the file:

    ```
    export PYTHONPATH=$PYTHONPATH:/home/groups/rhiju/rkretsch/PK
    export ARNIEFILE="/home/groups/rhiju/rkretsch/PK/arnie/arnie_file.txt"
    export NUPACKHOME="/home/groups/rhiju/rkretsch/PK/nupack3.0.6"
    ```
    - close and save the file by hitting control+x, then y to save the file, then enter to exit

5. on sherlock, navigate to your folder in scratch and create a directory in which to run viralknots

    ```
    cd GROUP_SCRATCH
    cd <sunetid>
    mkdir viralknots_run1
    cd viralknots_run1
    ```
    - now run the command `pwd` to get the path to the folder (you will use this path to write your sbatch files and copy your fasta file over to the cluster)

6. Now you must copy your fasta file over to sherlock to run in the folder you have just created
    - in the terminal on your home computer, navigate to the folder where you have saved your fasta file and run the following command:

    ```
    scp <fasta_file> <sunetid>@login.sherlock.stanford.edu:<path_to_viralknots_run1_folder>
    ```

7. Now you are ready to create your sbatch files to run your job on Viralknots. There are two sbatch files you must write: if you wish to run a large number of child jobs in parallel (which will greatly speed up compute time), you must create a template sbatch file for those jobs; you must also create a general sbatch file to run the parent job
    - first let's write the template sbatch file. Create the file as follows:

    `nano template_sbatch`

    - this command will take you to an empty file, in which you should write the following:

    
    ```
        #!/bin/bash

        #SBATCH -J viralknots_single #name of the job
        #SBATCH -o /<path_to_viralknots_run1>/single.out #file for job to write outputs
        #SBATCH -e /<path_to_viralknots_run1>/single.err #file for job to write errors
        #SBATCH -p biochem #partition to run on
        #SBATCH -t 7:00:00 #time for job to run
        #SBATCH -n 1 #number of cpus
        #SBATCH -N 1 #number of nodes

        cd /<path_to_viralknots_run1>
        module load gcc #if you want to run rnastructure
        module load glpk #if you want to run ipknots
        module load mpfr #if you want to run ipknots
    ```

    - Remember, DO NOT actually write a line with a command here. This will be done by the pipeline.
    - Now save the file by hitting control+x, then y to save, then enter to exit

8. Now you must create a sbatch file that you will use to submit the job with your desired arguments. Create the file using the following command:

    `nano viralknots_run1_sbatch`

    - write the file using the following template:

    ```
        #!/bin/bash

        #SBATCH -J viralknots_run1 #Name of job
        #SBATCH -o /<path_to_viralknots_run1/viralknots_run1.out #file where output is written
        #SBATCH -e /<path_to_viralknots_run1>/viralknots_run1.err #file where errors are written
        #SBATCH -p biochem #partition to run on (you should use biochem for now)
        #SBATCH -t 7:00:00 #time to run the job
        #SBATCH -n 1 #number of cpus
        #SBATCH -N 1 #number of gpus

        cd /<path_to_viralknots_run1> #go to wherever you want to run the job
        module load gcc #for rnastructure
        module load glpk #for ipknots
        module load mpfr #for ipknots

        python /home/groups/rhiju/rkretsch/PK/arnie/scripts/ViralKnots/ViralKnots.py <--pk_predict> <--shapeknots> -s <seq_filename> --step <step> -w <window> --pk_predictors <list_of_pk_predictors>  --bpp_package <bpp_package> --shape_data_folder <path_to_folder> --shape_data_sets <list_of_data_sets> <--shape_rankings> <--spawn> --template_sbatch <template_sbatch> --num_jobs <num_jobs>
    ```

    - the necessary inputs are as follows:
        - pk_predict: use this command to specify if you want to run pk predictors other than shapeknots; do not include this command if you only want to run shapeknots
        - shapeknots: use this command to specify if you want to run shapeknots; do not include this if you only want to run pk predictors
        - seq_filename: the name and location of the fasta file with the viral sequence
        - step: the number of nucleotides for the pipeline to slide down before starting the next window
        - window: the size of the nucleotide window to run predictions over
        - pk_predictors: a list of names of pk predictors to use; these should be formatted as a list separated by spaces NO COMMAS (e.g.: threshknot spotrna knotty pknots)
        - bpp_package: name of the bpp package you want to use if you are running threshknot; default is contrafold2
        - shape_data_folder: default is None; the location of the folder in which you have stored csv's with shape data
        - shape_data_sets: default is None; your shape data should be in files ending with .csv and the shape_data_sets variable is the names of the files WITHOUT the .csv
            - note that the format of the shape data in the file should be one reactivity value per line with no commas or any other comments/labels
        - shape_rankings: include this command if you want to receive scores for each structure based on agreement with available shape data
        - spawn: include this command if you want to spawn child jobs to run in parallel; recommended if you want to speed up processing time
        - template_sbatch: the name and location of the template sbatch file you created in step 5
        - num_jobs: the number of child jobs you want to spawn (max is ~1000 per hour on sherlock)
            - if running shapeknots, this needs to be at least as many as the number of shape data sets you are using

    - now save the file by typing control+x, then y to save changes, then enter to exit

6. submit the job by running the following command:
    
    `sbatch viralknots_sbatch`

- you can check the status of your job on sherlock by running the following:
    
    `squeue -u <sunetid>`

7. the output of your file will be a csv containing the name of the predictor used, the start and end location, the sequence, the predicted secondary structure in dotbracket notation, a True/False column denoting whether or not the structure represents a pseudoknot, and columns with the average F1 score and shape consensus scores for the entire structure as well as only the pseudoknotted base pairs
