-   **BaseBatchScripts** A generally functional version of the shell
    scripts used initiate the pyABC based parameter inference on an HPC.
    When executed, it starts a Redis-server, connects workers to the
    server and runs the python program containing the inference
    assignments. Originally provided by Emad Alamoudi.

-   **Models** Contains all files required to run the models, evaluate
    their results and create the figures included in the body of this
    thesis. Also contains the results of the runs we performed, however
    mostly without the database in which the full run details are saved,
    as these files are too large (usually several hundred MB for each
    run). The databases remain saved on the cluster infrastructure on
    which the respective test was run and can be made available upon
    request.

    -   **M1\_Tumor** Everything required to run and evaluate the tests
        of the tumor growth model (M1). The inference is performed in
        the python programs *TumorAdaptiveEps.py* and *TumorListEps.py*
        with an adaptive and a static epsilon respectively. Executing
        these files returns detailed information about each run in a
        database file, and some additional statistics about the effect
        of the look-ahead scheduling is written into .csv file.

        -   **BatchScript** A version of the shell scripts tailored to
            the tumor model test runs. The main changes are a different
            worker call, necessary to prevent daemonic processes, and
            the amount of active nodes being handed to the python
            program call as an argument.

        -   **Figures** The jupyter notebooks load the run data in order
            to visualize the results as seen in Figures 17-22 (created
            in *TumorVisualization.py*) and Figure 23 (created in
            *TumorVisStat.py*). The figures are the ones summarizing the
            results of the runs with a static epsilon schedule on equal
            conditions (see Section [sec:TumorRT]).

        -   **Testresults** The results of the performed (M1) runs,
            sorted by population size and workers used. Does not include
            the databases themselves, but instead some extracted
            information summarizing each generation, the .csv files and
            all visualizations for each individual run.

    -   **M2\_HIV** The python programs used to perform parameter
        inference for the HIV model (M2) and the results of the 8 runs
        (.csv files, and for the LA, adaptive epsilon run on 256 workers
        also the database). Also, the notebook used to visualize these
        results, creating Figures 25, 26. Model and posterior plots
        (Figure 24) provided by Nils Bundgaard.

    -   **T1\_ODE** Everything required to run and evaluate the tests of
        the ODE model (T1). The inference is performed multiple times by
        executing the python program *ODEWLogfiles.py*. Executing this
        file returns detailed information about each repetition in a
        separate database file, some additional statistics about the
        effect of the look-ahead scheduling for each run in different
        .csv files and additionally a .txt file containing summarized
        wall time statistics.

        -   **BatchScript** A version of the shell scripts tailored to
            the ODE model test runs. Mainly features some additional
            arguments being passed between the shell scripts.

        -   **Figures** The python program (*ODEBoxPlotCreation.py*)
            reads in all databases and .csv files to plot the
            summarizing graphs showing equivalence of the posteriors as
            seen in Figures 3-5. The jupyter notebook creates the wall
            time comparison plots (Figures 6-10) from the wall time
            summaries in the .txt files.

        -   **Testresults** The summarized wall time results of the test
            sorted by the runtime variance and the amount of nodes.
            Databases and .csv logfiles remain on the server, as there
            are several hundred of them.

    -   **T2\_MJP** Everything required to run and evaluate the tests of
        the Gillespie algorithm based model (T2). Basically a copy of
        the files in T1\_ODE, only with adapted paths and the model in
        the test file (*MJPRuntimeTest.py*) changed to the one for (M2).
        See above for details about the sub-folders.

    -   **T3\_UnbModes** Files used to run model (T3) both locally
        (*UnbalancedModes.ipynb*) and with minor modifications on an HPC
        (*UnbalancedModes.py*). Additionally the returned databases of
        the HPC runs (once using DYN, once LA scheduling) and the
        visualizations of the local and the cluster runs.

-   **Other Figures** The remaining Figures used in the main body
    (Figure 1,2 and the sketch of the tumor model). Including the
    program used to create Figure [fig:strategies], which reads in a
    file with a list of accepted or rejected particles with their
    simulation times and imitates the an STAT, DYN and LA scheduling
    approach to show how these particles would be distributed on the
    workers.

-   **pyABC** Two versions of the pyABC package: Version 0.10.14 which
    was used for testing and a modified version, including some changes
    to the built-in visualization functions and a first implementation
    of some of the Enhancements mentioned in Section [sec:OPTESS].

-   **tumor2d** Version 1.0.0 of the tumor2d package required to run the
    tumor model (M1).


