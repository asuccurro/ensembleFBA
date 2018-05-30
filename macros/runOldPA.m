tic
fprintf('Problem set up     ... success\n');
fprintf('Starting gap fill  ... (should finish in roughly 100 seconds)\n');
[modelList1] = build_network(old_seed_rxns_mat,biologicalData,params);
time2run = toc;
fprintf('Gap fill complete  ... success (%1.1f seconds)\n',time2run);
