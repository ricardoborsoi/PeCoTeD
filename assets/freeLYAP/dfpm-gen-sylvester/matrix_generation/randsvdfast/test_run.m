warning('off', 'randsvd_unit_tests:invalidUsage');
test_results = run(randsvdfast_unit_tests);
warning('on', 'randsvd_unit_tests:invalidUsage');
table(test_results)
