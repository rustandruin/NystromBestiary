% Runs all the experiments needed to generate the plots in the paper

% Apparently this special variable will be deprecated soon
maxNumCompThreads = 1; % run everything serially for accurate timing data
sendEmails = false; % send emails at the start and end of each experiment 

if (sendEmails)
    % note you need to setpref('Internet', 'SMTP_Password', 'pwd') manually
    % before calling this function!
    emailaddress = 'dummy@dummy.com';
    smtpserver = 'smtp.gmail.com';
    setpref('Internet', 'E_mail', emailaddress);
    setpref('Internet', 'SMTP_Username', emailaddress);
    setpref('Internet', 'SMTP_Server', smtpserver);

    % these settings taken from
    % http://www.mathworks.com/help/matlab/ref/sendmail.html
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth', 'true');
    props.setProperty('mail.smtp.socketFactory.class', ...
                        'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port', '465');

    notifier(emailaddress, @run_laplacians)
    notifier(emailaddress, @run_linear)
    notifier(emailaddress, @run_alg1)
    notifier(emailaddress, @run_rbf)
    notifier(emailaddress, @run_compact_rbf)
else
    run_laplacians
    run_linear
    run_alg1
    run_rbf
    run_compact_rbf
end
