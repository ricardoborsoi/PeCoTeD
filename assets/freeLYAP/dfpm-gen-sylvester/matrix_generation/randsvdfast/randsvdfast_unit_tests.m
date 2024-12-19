function tests = randsvdfast_unit_tests
%FUN_BASED_UNIT_TESTS   Function-based unit tests for randsvdfast.
%   This function is not meant to be called directly.
%   To run the tests, please use the script run_test instead.
  warning('randsvd_unit_tests:invalidUsage',...
          'Please use the script ''test_run'' to run the tests.');
  tests = functiontests(localfunctions);
end

function testDefaults(testCase)
  rng(1); A = randsvdfast();
  rng(1); B = randsvdfast(10,10,3,1,0,'double',true);
  testCase.verifyEqual(A,B);
end

function testInvalidN(testCase)
  testfun = @(x)randsvdfast(x);
  for n = [-Inf,-1.1,-1,0,1.1,Inf,NaN]
    testCase.verifyError(@()testfun(n), 'randsvd:invalidN');
  end
  testCase.verifyError(@()testfun({}), 'randsvd:invalidN');
  testCase.verifyError(@()testfun([]), 'randsvd:invalidN');
  testCase.verifyError(@()testfun([1,2,3]), 'randsvd:invalidN');
  testCase.verifyError(@()testfun('c'), 'randsvd:invalidN');
  testCase.verifyError(@()testfun('string'), 'randsvd:invalidN');
  testCase.verifyError(@()testfun(@(x)(x)), 'randsvd:invalidN');
end

function testInvalidKAPPA(testCase)
  n = 3;
  testfun = @(x)randsvdfast(n,x);
  for kappa = [-Inf,-1.1,-1,0,0.1,Inf,NaN]
    testCase.verifyError(@()testfun(kappa), 'randsvd:invalidKAPPA');
  end
  testCase.verifyError(@()testfun({}), 'randsvd:invalidKAPPA');
  testCase.verifyError(@()testfun([]), 'randsvd:invalidKAPPA');
  testCase.verifyError(@()testfun('c'), 'randsvd:invalidKAPPA');
  testCase.verifyError(@()testfun('string'), 'randsvd:invalidKAPPA');
  testCase.verifyError(@()testfun(@(x)(x)), 'randsvd:invalidKAPPA');
  testCase.verifyError(@()testfun(rand(1,n-1)), 'randsvd:invalidKAPPA');
  testCase.verifyError(@()testfun(rand(1,n+1)), 'randsvd:invalidKAPPA');

  testCase.verifyError(@()randsvdfast(n,10,6), 'randsvd:invalidKappaMode');
end

function testInvalidMODE(testCase)
  testfun = @(x)randsvdfast(2,10,x);
  for mode = [-Inf,-1.1,-1,0.1,7,Inf,NaN]
    testCase.verifyError(@()testfun(mode), 'randsvd:invalidMODE');
  end
  testCase.verifyError(@()testfun({}), 'randsvd:invalidMODE');
  testCase.verifyError(@()testfun([]), 'randsvd:invalidMODE');
  testCase.verifyError(@()testfun([1,2]), 'randsvd:invalidMODE');
  testCase.verifyError(@()testfun('c'), 'randsvd:invalidMODE');
  testCase.verifyError(@()testfun('string'), 'randsvd:invalidMODE');
  testCase.verifyError(@()testfun(@(x)(x)), 'randsvd:invalidMODE');
end

function testInvalidMETHOD(testCase)
  n = 2;
  testfun = @(x)randsvdfast(n,10,3,x);
  for method = [-Inf,-1.1,-1,0.1,5,Inf,NaN]
    testCase.verifyError(@()testfun(method), 'randsvd:invalidMETHOD');
  end
  testCase.verifyError(@()testfun({}), 'randsvd:invalidMETHOD');
  testCase.verifyError(@()testfun([]), 'randsvd:invalidMETHOD');
  testCase.verifyError(@()testfun([1,2]), 'randsvd:invalidMETHOD');
  testCase.verifyError(@()testfun('c'), 'randsvd:invalidMETHOD');
  testCase.verifyError(@()testfun('string'), 'randsvd:invalidMETHOD');
  testCase.verifyError(@()testfun(@(x)(x)), 'randsvd:invalidMETHOD');
  testCase.verifyError(@()randsvdfast(n,[1:n],3,3), 'randsvd:invalidKappaMethod');
  testCase.verifyError(@()randsvdfast(n,[1:n],3,4), 'randsvd:invalidKappaMethod');
  for mode = [3:5]
    for method = [3:4]
      testCase.verifyError(@()randsvdfast(n,10,mode,method), 'randsvd:invalidModeMethod');
    end
  end
  for mode = [0:2]
    for method = [3:4]
      testCase.verifyError(@()randsvdfast([2,3],10,mode,method),'randsvd:invalidNMethod');
      testCase.verifyError(@()randsvdfast([2,3],10,mode,method),'randsvd:invalidNMethod');
    end
  end
end

function testInvalidMATRIX(testCase)
  testfun = @(x)randsvdfast(2,10,3,1,x);
  for matrix = [-Inf,-1.1,-1,0.1,8,Inf,NaN]
    testCase.verifyError(@()testfun(matrix), 'randsvd:invalidMATRIX');
  end
  testCase.verifyError(@()testfun({}), 'randsvd:invalidMATRIX');
  testCase.verifyError(@()testfun([]), 'randsvd:invalidMATRIX');
  testCase.verifyError(@()testfun([1,2]), 'randsvd:invalidMATRIX');
  testCase.verifyError(@()testfun('c'), 'randsvd:invalidMATRIX');
  testCase.verifyError(@()testfun('string'), 'randsvd:invalidMATRIX');
end

function testInvalidCLASSNAME(testCase)
  testfun = @(x)randsvdfast(2,10,3,1,0,x);
  for classname = [-Inf,-1,0,0.1,1,Inf,NaN]
    testCase.verifyError(@()testfun(classname), 'randsvd:invalidCLASSNAME');
  end
  testCase.verifyError(@()testfun({}), 'randsvd:invalidCLASSNAME');
  testCase.verifyError(@()testfun([]), 'randsvd:invalidCLASSNAME');
  testCase.verifyError(@()testfun([1,2]), 'randsvd:invalidCLASSNAME');
  testCase.verifyError(@()testfun('c'), 'randsvd:invalidCLASSNAME');
  testCase.verifyError(@()testfun('mp'), 'randsvd:invalidCLASSNAME');
  testCase.verifyError(@()testfun('string'), 'randsvd:invalidCLASSNAME');
end

function testInvalidREALOUT(testCase)
  testfun = @(x)randsvdfast(2,10,3,1,0,'double',x);
  for realout = [-Inf,-1,0,0.1,1,Inf,NaN]
    testCase.verifyError(@()testfun(realout), 'randsvd:invalidREALOUT');
  end
  testCase.verifyError(@()testfun({}), 'randsvd:invalidREALOUT');
  testCase.verifyError(@()testfun([]), 'randsvd:invalidREALOUT');
  testCase.verifyError(@()testfun([1,2]), 'randsvd:invalidREALOUT');
  testCase.verifyError(@()testfun('c'), 'randsvd:invalidREALOUT');
  testCase.verifyError(@()testfun('mp'), 'randsvd:invalidREALOUT');
  testCase.verifyError(@()testfun('string'), 'randsvd:invalidREALOUT');
end

function testInvalidArgumentNumber(testCase)
  testCase.verifyError(@()randsvdfast(2,10,3,1,0,'double',true,1), 'MATLAB:TooManyInputs');
end

function testSingularValues(testCase)
  classnames = {'single', 'double'};
  for i = 1:2
    classname = classnames{i};
    if strcmp(classname, 'single');
      exponents = [1 2 3 4 5];
    else
      exponents = [1 2 4 8 12];
    end
    for n = 2.^[2 3 4]
      for kappa = 10.^exponents
        for mode = 0:6
          if mode < 3
            maxmethod = 4;
          else
            maxmethod = 2;
          end
          rng(1);
          switch mode
            case 0
              sigma = [1; ones(n-2, 1)./sqrt(kappa); 1./(kappa)];
            case 1
              sigma = [kappa; ones(n-1, 1)]./kappa;
            case 2
              sigma = [ones(n-1, 1); 1/kappa];
            case 3
              sigma = nthroot(kappa, 1-n).^(0:n-1)';
            case 4
              h = (kappa-1)/(n-1);
              sigma = fliplr(1 + h*(0:n-1))' / kappa;
            case 5
              sigma = sort(exp(-rand(n, 1)*log(kappa)), 'descend');
            case 6
              kappa = [n:-1:1]';
              sigma = [n:-1:1]';
          end
          for method = 1:maxmethod
            if method > 2
              if mode == 0
                sigmaloc = sigma * sqrt(kappa);
              elseif mode == 1
                sigmaloc = sigma * kappa;
              end
            else
              sigmaloc = sigma;
            end
            for matrix = 0:7
              for realout = [false, true]
                rng(1);
                A = randsvdfast(n,kappa,mode,method,matrix,classname,realout);
                if norm(svd(A)-sigmaloc)./norm(sigmaloc) > sqrt(eps(classname))
                  [svd(A), sigmaloc]
                  keyboard
                end
                testCase.verifySize(A,[n,n]);
                % testCase.verifyLessThan(abs(cond(A)-kappa)/kappa,0.02);
                testCase.verifyLessThan(norm(svd(A)-sigmaloc)./norm(sigmaloc),0.02);

                if method < 3
                  rng(1);
                  A = randsvdfast([2*n,n],kappa,mode,method,matrix,classname,realout);
                  testCase.verifySize(A,[2*n,n]);
                  testCase.verifyLessThan(norm(svd(A)-sigmaloc)./norm(sigmaloc),0.02);
                  % testCase.verifyLessThan(abs(cond(A)-kappa)/kappa,0.02);

                  rng(1);
                  A = randsvdfast([n,2*n],kappa,mode,method,matrix,classname,realout);
                  testCase.verifySize(A,[n,2*n]);
                  testCase.verifyLessThan(norm(svd(A)-sigmaloc)./norm(sigmaloc),0.02);
                  % testCase.verifyLessThan(abs(cond(A)-kappa)/kappa,0.02);
                end
              end
            end
          end
        end
      end
    end
  end
end
