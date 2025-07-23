function mailme(title,message,file_pos)
mail = 'matlab_nzj@163.com';  % ①邮箱地址
password = 'AFesDWm5QEcRfBzn'; % ②密码

% 服务器设置
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.163.com'); % ③SMTP服务器
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');



% 收件人
receiver='nizijun@stu.pku.edu.cn'; 
% 邮件标题
mailtitle=title;
% 邮件内容
mailcontent=message;
% 发送
if(file_pos == 0)
    sendmail(receiver, mailtitle, mailcontent);
else
    sendmail(receiver, mailtitle, mailcontent, file_pos);
end
end