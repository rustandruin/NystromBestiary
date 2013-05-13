% Runs all the experiments needed to generate the plots in the paper

% Apparently this special variable will be deprecated soon
maxNumCompThreads = 1; % run everything serially for accurate timing data

% Send notifications to your email address upon start and 
% completion of tasks
% uncomment if you have sendmail configured with matlab
%
%emailaddress = 'dummy@dummy.com';
%
%notifier(emailaddress, @run_laplacians)
%notifier(emailaddress, @run_linear)
%notifier(emailaddress, @run_alg1)
%notifier(emailaddress, @run_rbf)
%notifier(emailaddress, @run_compact_rbf)

run_laplacians
run_linear
run_alg1
run_rbf
run_compact_rbf

