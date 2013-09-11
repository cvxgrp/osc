% Creates a local directory called cvxgen_CODENUM, and places a solver called
% csolve in your local directory.

% TODO: make platform independent.

function cvxgen(code)

code = num2str(code);

url = ['http://cvxgen.stanford.edu/matlab_test/' code];
if ~strcmp(urlread(url), 'success')
  disp(['Failed to retrieve problem ' code '.']);
  return;
else
  disp('Retrieving solver from http://cvxgen.stanford.edu/...');
end

url = ['http://cvxgen.stanford.edu/matlab/' code];
dir = 'cvxgen/';

if exist(dir) ~= 7
  status = mkdir(dir);
end

unzip(url, '.');

if exist([dir '/solver.c']) == 2
  disp(['Downloaded to ' dir '.']);
else
  disp('Failed to retrieve file. Giving up.');
  return;
end

disp(' ');
disp('cvxgen is in beta. Please note that there is no warranty.');
disp('You may want to check csolve(params) against cvxsolve(params), which uses cvx.');
disp(' ');
disp('Compiling...');

cd(dir);
try
  make_csolve;
  success = 1;
catch exception
  disp('!!! Compilation failed.')
  cd('..');
  success = 0;
  rethrow(exception);
end

if success
  pause(0.05); % Attempt to avoid a certain race condition.
  cd('..');
  pause(0.05);

  pause(0.05); % Attempt to avoid a certain race condition.
  copyfile([dir '/csolve.m*'], '.');
  copyfile([dir '/cvxsolve.m'], '.');

  disp('Success. Type help csolve, or')
  disp(' ');
  disp('  [vars, status] = csolve(params, settings)');
  disp(' ');
end


