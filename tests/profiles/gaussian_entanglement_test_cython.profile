         225479 function calls (223838 primitive calls) in 1.639 seconds

   Ordered by: internal time

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
    50288    1.456    0.000    1.456    0.000 {pyge.cythonlib.cython_gaussian_entanglement.cython_gaussian_entanglement}
        1    0.078    0.078    1.549    1.549 /home/leonardo/pyge/pyge/gent.py:79(ge_loops)
        1    0.018    0.018    0.038    0.038 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/PDBParser.py:224(_parseatoms)
   100890    0.016    0.000    0.016    0.000 {built-in method builtins.abs}
        1    0.005    0.005    0.006    0.006 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/PDB.py:377(_read_frame)
        9    0.004    0.000    0.005    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:683(__init__)
    212/2    0.004    0.000    0.010    0.005 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:494(_parse)
    22462    0.003    0.000    0.003    0.000 {method 'append' of 'list' objects}
    12853    0.003    0.000    0.003    0.000 {method 'strip' of 'str' objects}
        2    0.003    0.001    0.003    0.002 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/base.py:145(change_squash)
      157    0.002    0.000    0.006    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/numeric.py:2354(within_tol)
  674/176    0.002    0.000    0.018    0.000 {built-in method numpy.core._multiarray_umath.implement_array_function}
        1    0.002    0.002    0.009    0.009 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/PDB.py:249(__init__)
      157    0.002    0.000    0.016    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/numeric.py:2273(isclose)
      319    0.002    0.000    0.002    0.000 {method 'reduce' of 'numpy.ufunc' objects}
    266/2    0.002    0.000    0.003    0.002 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:87(_compile)
       25    0.002    0.000    0.002    0.000 {built-in method numpy.array}
     2520    0.002    0.000    0.002    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/PDBParser.py:90(float_or_default)
      314    0.001    0.000    0.003    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/_ufunc_config.py:32(seterr)
      315    0.001    0.000    0.005    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/fromnumeric.py:69(_wrapreduction)
      314    0.001    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/_ufunc_config.py:131(geterr)
     3063    0.001    0.000    0.001    0.000 {method 'startswith' of 'str' objects}
   287/23    0.001    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:175(getwidth)
        1    0.001    0.001    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/PDBParser.py:429(<genexpr>)
      312    0.001    0.000    0.003    0.000 {method 'all' of 'numpy.generic' objects}
     1421    0.001    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:165(__getitem__)
      903    0.001    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:234(__next)
       90    0.001    0.000    0.002    0.000 <frozen importlib._bootstrap_external>:1536(find_spec)
      315    0.001    0.000    0.005    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/fromnumeric.py:2406(all)
        1    0.001    0.001    1.639    1.639 /home/leonardo/pyge/tests/gent_test.py:13(gaussian_entanglement_test)
      476    0.001    0.000    0.001    0.000 {built-in method numpy.asanyarray}
     47/2    0.001    0.000    0.010    0.005 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:436(_parse_sub)
2536/2224    0.001    0.000    0.001    0.000 {built-in method builtins.len}
     1261    0.001    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/PDBParser.py:425(<genexpr>)
        1    0.001    0.001    0.003    0.003 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/PDBParser.py:410(_parsebonds)
      315    0.000    0.000    0.007    0.000 <__array_function__ internals>:177(all)
       91    0.000    0.000    0.000    0.000 {built-in method posix.stat}
        1    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/guessers.py:120(<listcomp>)
       88    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:239(_add_prop)
      630    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:255(get)
     2520    0.000    0.000    0.000    0.000 {method 'capitalize' of 'str' objects}
        1    0.000    0.000    0.041    0.041 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/PDBParser.py:200(parse)
      314    0.000    0.000    0.000    0.000 {built-in method numpy.seterrobj}
     1260    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/guessers.py:441(get_atom_mass)
     1666    0.000    0.000    0.000    0.000 {built-in method builtins.isinstance}
      450    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap_external>:128(<listcomp>)
      450    0.000    0.000    0.001    0.000 <frozen importlib._bootstrap_external>:126(_path_join)
     1517    0.000    0.000    0.000    0.000 {method 'readline' of '_io.BufferedReader' objects}
      157    0.000    0.000    0.002    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/_ufunc_config.py:429(__enter__)
      157    0.000    0.000    0.002    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/_ufunc_config.py:434(__exit__)
      315    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/fromnumeric.py:70(<dictcomp>)
      314    0.000    0.000    0.002    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/_methods.py:60(_all)
       51    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:292(_optimize_charset)
       19    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:892(_process_attr)
        1    0.000    0.000    0.000    0.000 {method 'splitlines' of 'str' objects}
      486    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:250(match)
      157    0.000    0.000    0.017    0.000 <__array_function__ internals>:177(isclose)
      289    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:112(__init__)
      628    0.000    0.000    0.000    0.000 {built-in method numpy.geterrobj}
       18    0.000    0.000    0.003    0.000 <frozen importlib._bootstrap_external>:1399(_get_spec)
      409    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:173(append)
      721    0.000    0.000    0.000    0.000 {built-in method builtins.min}
      157    0.000    0.000    0.001    0.000 <__array_function__ internals>:177(result_type)
        8    0.000    0.000    0.000    0.000 {built-in method io.open}
      446    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:161(__len__)
        1    0.000    0.000    0.000    0.000 {built-in method _pickle.load}
        7    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/base.py:182(<genexpr>)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:210(_mix)
      214    0.000    0.000    0.000    0.000 {method 'decode' of 'bytes' objects}
       18    0.000    0.000    0.003    0.000 <frozen importlib._bootstrap>:921(_find_spec)
       42    0.000    0.000    0.013    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/tokenize.py:431(_tokenize)
      157    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/_ufunc_config.py:425(__init__)
       33    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/codecs.py:319(decode)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/base.py:179(get_borders)
        1    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:606(_apply)
       48    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:356(_escape)
    28/18    0.000    0.000    0.004    0.000 <frozen importlib._bootstrap>:1022(_find_and_load)
       15    0.000    0.000    0.000    0.000 {method 'read' of '_io.BufferedReader' objects}
      336    0.000    0.000    0.000    0.000 {built-in method builtins.getattr}
       79    0.000    0.000    0.001    0.000 {built-in method builtins.any}
      450    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:244(_verbose_message)
       56    0.000    0.000    0.000    0.000 {method 'format' of 'str' objects}
      112    0.000    0.000    0.000    0.000 {method 'match' of 're.Pattern' objects}
        8    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/picklable_file_io.py:196(__init__)
      199    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:287(tell)
      324    0.000    0.000    0.000    0.000 {method 'items' of 'dict' objects}
       33    0.000    0.000    0.000    0.000 {built-in method _codecs.utf_8_decode}
      315    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/fromnumeric.py:2401(_all_dispatcher)
      900    0.000    0.000    0.000    0.000 {method 'rstrip' of 'str' objects}
        2    0.000    0.000    0.013    0.007 /home/leonardo/.conda/envs/ge/lib/python3.10/re.py:288(_compile)
       28    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:179(_get_module_lock)
       50    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:433(_uniq)
       28    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:125(release)
        2    0.000    0.000    0.000    0.000 {built-in method builtins.compile}
        8    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:176(_subclass)
        1    0.000    0.000    0.000    0.000 {built-in method numpy.fromfile}
       11    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/numeric.py:76(zeros_like)
      459    0.000    0.000    0.000    0.000 {method 'join' of 'str' objects}
    28/18    0.000    0.000    0.004    0.000 <frozen importlib._bootstrap>:987(_find_and_load_unlocked)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/arraysetops.py:523(in1d)
        1    0.000    0.000    0.000    0.000 {method 'outer' of 'numpy.ufunc' objects}
      113    0.000    0.000    0.000    0.000 {built-in method builtins.setattr}
      373    0.000    0.000    0.000    0.000 {built-in method builtins.ord}
        9    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:429(_get_stream)
       15    0.000    0.000    0.000    0.000 {built-in method numpy.asarray}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/tokenize.py:185(untokenize)
        1    0.000    0.000    0.000    0.000 {MDAnalysis.lib._cutil.unique_int_1d}
       51    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:265(_compile_charset)
      148    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:614(<genexpr>)
        1    0.000    0.000    0.002    0.002 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:153(_generate_from_topology)
        2    0.000    0.000    0.013    0.007 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/format.py:547(_filter_header)
       28    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:100(acquire)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/picklable_file_io.py:108(__init__)
        1    0.000    0.000    0.000    0.000 {method 'sort' of 'numpy.ndarray' objects}
       74    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/fnmatch.py:70(fnmatchcase)
      305    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/PDBParser.py:394(<genexpr>)
       19    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topology.py:510(add_TopologyAttr)
      157    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/numeric.py:2269(_isclose_dispatcher)
       50    0.000    0.000    0.000    0.000 {built-in method fromkeys}
        2    0.000    0.000    0.014    0.007 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/format.py:582(_read_array_header)
        8    0.000    0.000    0.000    0.000 {method 'close' of '_io.TextIOWrapper' objects}
       28    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:71(__init__)
      165    0.000    0.000    0.000    0.000 {method 'get' of 'dict' objects}
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_compression.py:39(__init__)
       18    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/_distutils_hack/__init__.py:89(find_spec)
      108    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap_external>:1356(_path_importer_cache)
        2    0.000    0.000    0.014    0.007 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/npyio.py:263(load)
      185    0.000    0.000    0.000    0.000 {built-in method builtins.max}
       21    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:85(opengroup)
      191    0.000    0.000    0.000    0.000 {method 'find' of 'bytearray' objects}
       12    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:279(__setattr__)
       13    0.000    0.000    0.000    0.000 {method 'close' of '_io.BufferedReader' objects}
      211    0.000    0.000    0.000    0.000 {method 'strip' of 'bytes' objects}
       90    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap_external>:140(_path_stat)
        1    0.000    0.000    0.056    0.056 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:317(__init__)
       11    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(zeros_like)
        5    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:56(parse_parts)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:576(__getitem__)
        3    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:324(anyopen)
       57    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:447(_simple)
       19    0.000    0.000    0.000    0.000 {built-in method numpy.zeros}
       12    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(copyto)
        1    0.000    0.000    0.000    0.000 /home/leonardo/pyge/pyge/gent.py:40(_midposition_vectors)
       28    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:198(cb)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topology.py:451(__init__)
      7/5    0.000    0.000    0.000    0.000 {built-in method _abc._abc_subclasscheck}
        1    0.000    0.000    0.010    0.010 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:488(load_new)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/arraysetops.py:323(_unique1d)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/bz2.py:37(__init__)
      157    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/multiarray.py:664(result_type)
        6    0.000    0.000    0.000    0.000 {method 'readline' of '_io.TextIOWrapper' objects}
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:134(__init__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/pyge/pyge/gent.py:17(_loop_list)
        2    0.000    0.000    0.013    0.006 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:783(compile)
       58    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:169(__setitem__)
       21    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:97(closegroup)
       28    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:169(__enter__)
        2    0.000    0.000    0.003    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/_get_readers.py:31(get_reader_for)
        5    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:569(_parse_args)
        5    0.000    0.000    0.000    0.000 {built-in method numpy.arange}
     18/2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/ast.py:84(_convert)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:128(make_classes)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:533(__exit__)
       11    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(empty_like)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:2744(positions)
       72    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:897(__exit__)
      104    0.000    0.000    0.000    0.000 {built-in method builtins.hasattr}
        2    0.000    0.000    0.014    0.007 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/format.py:719(read_array)
        7    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:405(__init__)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_compression.py:66(readinto)
        2    0.000    0.000    0.000    0.000 {built-in method _sre.compile}
       46    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:82(groups)
       72    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:893(__enter__)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:560(_compile_info)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:226(int_array_is_sorted)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:435(_mk_bitmap)
       38    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/tokenize.py:172(add_whitespace)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:437(<listcomp>)
        1    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:2862(select_atoms)
       18    0.000    0.000    0.003    0.000 <frozen importlib._bootstrap_external>:1431(find_spec)
      118    0.000    0.000    0.000    0.000 {method 'rpartition' of 'str' objects}
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:547(__init__)
       20    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/logging/__init__.py:1724(isEnabledFor)
       88    0.000    0.000    0.000    0.000 {method 'add' of 'set' objects}
        7    0.000    0.000    0.000    0.000 {built-in method _abc._abc_instancecheck}
        5    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:621(__str__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/posixpath.py:400(_joinrealpath)
        2    0.000    0.000    0.000    0.000 {built-in method _functools.reduce}
      128    0.000    0.000    0.000    0.000 {built-in method _imp.acquire_lock}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/ast.py:54(literal_eval)
      128    0.000    0.000    0.000    0.000 {built-in method _imp.release_lock}
        4    0.000    0.000    0.000    0.000 {built-in method posix.lstat}
      2/1    0.000    0.000    0.009    0.009 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:2479(wrapper)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_compression.py:59(close)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:1569(wrapper)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/inspect.py:3046(_bind)
       28    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:173(__exit__)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/picklable_file_io.py:341(pickle_open)
       18    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:357(__init__)
        4    0.000    0.000    0.000    0.000 {method 'read1' of '_io.BufferedReader' objects}
        1    0.000    0.000    0.002    0.002 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/_get_readers.py:205(get_parser_for)
        2    0.000    0.000    0.010    0.005 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:944(parse)
     10/8    0.000    0.000    0.002    0.000 {built-in method builtins.__import__}
       20    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/logging/__init__.py:1455(debug)
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/posixpath.py:71(join)
        4    0.000    0.000    0.000    0.000 {method 'copy' of 'numpy.ndarray' objects}
       18    0.000    0.000    0.000    0.000 {built-in method _imp.is_builtin}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:1431(parse)
       23    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:81(_combine_flags)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_compression.py:72(read)
       18    0.000    0.000    0.000    0.000 {built-in method builtins.locals}
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/picklable_file_io.py:413(bz2_pickle_open)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/enum.py:986(__and__)
       90    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap_external>:71(_relax_case)
        2    0.000    0.000    0.000    0.000 {method 'nonzero' of 'numpy.ndarray' objects}
       49    0.000    0.000    0.000    0.000 {built-in method __new__ of type object at 0x5612ffd4e5c0}
       18    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:746(find_spec)
        8    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/codecs.py:309(__init__)
       40    0.000    0.000    0.000    0.000 <string>:1(<lambda>)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/picklable_file_io.py:264(__init__)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:615(_make_child)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/format.py:932(_read_bytes)
       56    0.000    0.000    0.000    0.000 {built-in method _thread.allocate_lock}
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:403(__init__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/twodim_base.py:376(tri)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/fnmatch.py:80(translate)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:791(_asunique)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:139(__exit__)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:430(_read_gzip_header)
        1    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/guessers.py:107(guess_masses)
       28    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:165(__init__)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:917(get_ext)
        2    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/converters/ParmEd.py:97(_format_hint)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topology.py:180(__init__)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:279(helper)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/genericpath.py:121(_splitext)
        1    0.000    0.000    0.043    0.043 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:104(_topology_from_file_like)
        5    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:608(_format_parsed_parts)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/twodim_base.py:497(triu)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/guessers.py:124(validate_atom_types)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/base.py:356(convert_pos_from_native)
        7    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:495(isstream)
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:698(universe)
       16    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_compression.py:12(_check_not_closed)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/picklable_file_io.py:485(gzip_pickle_open)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/bz2.py:97(close)
       56    0.000    0.000    0.000    0.000 {method '__exit__' of '_thread.lock' objects}
       18    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:826(find_spec)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:1505(_parse_subexp)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/utils.py:967(safe_eval)
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/enum.py:359(__call__)
        2    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/converters/OpenMM.py:83(_format_hint)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:589(_from_parts)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/posixpath.py:337(normpath)
       56    0.000    0.000    0.000    0.000 {built-in method _thread.get_ident}
        1    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:651(select_atoms)
        7    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:490(hasmethod)
        3    0.000    0.000    0.000    0.000 {built-in method zlib.decompressobj}
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:998(guess_format)
        1    0.000    0.000    0.000    0.000 {method 'argmin' of 'numpy.ndarray' objects}
        5    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:239(splitroot)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:482(enter_context)
       15    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:449(__len__)
       10    0.000    0.000    0.001    0.000 {built-in method builtins.next}
      2/1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:485(_get_literal_prefix)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:760(isunique)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:2446(wrapper)
        2    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/converters/OpenMM.py:154(_format_hint)
       12    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/multiarray.py:1071(copyto)
        4    0.000    0.000    0.001    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:273(openany)
        3    0.000    0.000    0.000    0.000 {method 'decompress' of '_bz2.BZ2Decompressor' objects}
       28    0.000    0.000    0.000    0.000 {built-in method posix.fspath}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/base.py:201(<listcomp>)
        5    0.000    0.000    0.000    0.000 {built-in method builtins.all}
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/picklable_file_io.py:326(__init__)
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:573(__len__)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:102(__init__)
     10/8    0.000    0.000    0.002    0.000 <frozen importlib._bootstrap>:233(_call_with_frames_removed)
        2    0.000    0.000    0.003    0.002 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:622(_code)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/getlimits.py:672(max)
        1    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:240(apply)
       11    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/base.py:187(<genexpr>)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/posixpath.py:117(splitext)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/converters/RDKit.py:139(_format_hint)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/converters/ParmEdParser.py:141(_format_hint)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:2494(__init__)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:86(read)
        7    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/abc.py:117(__instancecheck__)
        9    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:1045(iterable)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/fromnumeric.py:51(_wrapfunc)
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/enum.py:678(__new__)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:463(read)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:514(_push_exit_callback)
       12    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:296(_class_escape)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:225(__init__)
       11    0.000    0.000    0.000    0.000 {method 'extend' of 'list' objects}
       15    0.000    0.000    0.000    0.000 {built-in method sys.intern}
       20    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topology.py:498(n_atoms)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:3216(<genexpr>)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/format.py:286(descr_to_dtype)
        1    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/converters/OpenMMParser.py:232(_format_hint)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:755(issorted)
        2    0.000    0.000    0.000    0.000 {built-in method numpy.empty}
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:327(close)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/shape_base.py:285(hstack)
        2    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(nonzero)
       11    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/numeric.py:72(_zeros_like_dispatcher)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:306(__new__)
        3    0.000    0.000    0.000    0.000 {method 'astype' of 'numpy.ndarray' objects}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/picklable_file_io.py:150(__init__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:2702(asunique)
       22    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/bz2.py:121(closed)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/format.py:221(read_magic)
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/chain.py:370(_format_hint)
        8    0.000    0.000    0.000    0.000 {method 'replace' of 'str' objects}
        2    0.000    0.000    0.000    0.000 {method 'ravel' of 'numpy.ndarray' objects}
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:410(_init_read)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:971(format_from_filename_extension)
        1    0.000    0.000    0.000    0.000 /home/leonardo/pyge/pyge/gent.py:44(_bond_vectors)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/twodim_base.py:33(_min_int)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/converters/OpenMMParser.py:83(_format_hint)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/shape_base.py:23(atleast_1d)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:598(__init__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/fnmatch.py:44(_compile_pattern)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/tokenize.py:259(untokenize)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/arraysetops.py:138(unique)
       18    0.000    0.000    0.000    0.000 {built-in method _imp.is_frozen}
       38    0.000    0.000    0.000    0.000 {method 'span' of 're.Match' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/base.py:1352(__init__)
       10    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:705(ix)
        1    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(in1d)
        2    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(where)
        2    0.000    0.000    0.000    0.000 {method 'tell' of '_io.BufferedReader' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/inspect.py:2836(apply_defaults)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/logging/__init__.py:219(_acquireLock)
        2    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:130(__enter__)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/bz2.py:126(fileno)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/ast.py:33(parse)
        3    0.000    0.000    0.000    0.000 {method 'seek' of '_io.BufferedReader' objects}
        7    0.000    0.000    0.000    0.000 {method 'split' of 'str' objects}
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:937(check_compressed_format)
       20    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topology.py:502(n_residues)
       11    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/multiarray.py:80(empty_like)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/converters/RDKitParser.py:162(_format_hint)
        4    0.000    0.000    0.000    0.000 {method 'lstrip' of 'str' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:82(grab_not_keywords)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:928(fix_flags)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/bz2.py:166(read1)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/re.py:269(escape)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:868(_parse_flags)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:80(__init__)
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:631(__fspath__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:2520(_bondDict)
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:600(_from_parsed_parts)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:509(_push_cm_exit)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/units.py:346(get_conversion_factor)
        8    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/codecs.py:260(__init__)
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/os.py:1080(__subclasshook__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:456(_generate_overlap_table)
       20    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topology.py:506(n_segments)
       24    0.000    0.000    0.000    0.000 {method 'isidentifier' of 'str' objects}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/posixpath.py:60(isabs)
      7/5    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/abc.py:121(__subclasscheck__)
        7    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/bz2.py:135(readable)
        2    0.000    0.000    0.000    0.000 {built-in method _struct.calcsize}
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:853(__truediv__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/bz2.py:186(readline)
        6    0.000    0.000    0.000    0.000 {method 'cast' of 'memoryview' objects}
        6    0.000    0.000    0.000    0.000 {method 'translate' of 'bytearray' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/numeric.py:149(ones)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:303(read1)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/posixpath.py:391(realpath)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/posixpath.py:41(_get_sep)
        1    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(unique)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_parse.py:76(__init__)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/getlimits.py:659(min)
        1    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:217(_apply)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:1064(resolve)
        1    0.000    0.000    0.000    0.000 {method 'flatten' of 'numpy.ndarray' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:770(trajectory)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/bz2.py:131(seekable)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/inspect.py:191(isclass)
        4    0.000    0.000    0.000    0.000 {method 'readline' of '_io.StringIO' objects}
        4    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:619(isstring)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/fromnumeric.py:1866(nonzero)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/numeric.py:2194(allclose)
        6    0.000    0.000    0.000    0.000 {function DecompressReader.close at 0x7f157f076440}
       16    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:323(closed)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_collections_abc.py:78(_check_methods)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:442(_create_exit_wrapper)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:67(is_keyword)
        6    0.000    0.000    0.000    0.000 {method 'rfind' of 'str' objects}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:452(__init__)
        2    0.000    0.000    0.013    0.007 /home/leonardo/.conda/envs/ge/lib/python3.10/re.py:249(compile)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/H5MD.py:567(_format_hint)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/logging/__init__.py:1710(getEffectiveLevel)
       18    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/_distutils_hack/__init__.py:96(<lambda>)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/tokenize.py:166(__init__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:865(parent)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/format.py:633(<genexpr>)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/logging/__init__.py:228(_releaseLock)
        1    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(triu)
        2    0.000    0.000    0.000    0.000 {method 'translate' of 'str' objects}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/io.py:60(__getattr__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:2517(__len__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/MMTFParser.py:156(_format_hint)
        1    0.000    0.000    0.001    0.001 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:141(_resolve_coordinates)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:94(join_parsed_parts)
        1    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(concatenate)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/MMTF.py:62(_format_hint)
        4    0.000    0.000    0.000    0.000 {method 'startswith' of 'bytes' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:957(__new__)
        1    0.000    0.000    0.013    0.013 /home/leonardo/.conda/envs/ge/lib/python3.10/tokenize.py:99(_compile)
        2    0.000    0.000    0.000    0.000 {built-in method sys.exc_info}
      4/2    0.000    0.000    0.000    0.000 {method 'seekable' of '_io.BufferedReader' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/pathlib.py:1092(stat)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:759(trajectory)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_compression.py:16(_check_can_read)
        5    0.000    0.000    0.000    0.000 {method 'insert' of 'list' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/base.py:608(__exit__)
        3    0.000    0.000    0.000    0.000 {built-in method zlib.crc32}
        1    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(hstack)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/PDB.py:331(<listcomp>)
        1    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(allclose)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/posixpath.py:376(abspath)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/types.py:176(__get__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/shape_base.py:207(_arrays_for_stack_dispatcher)
        4    0.000    0.000    0.000    0.000 {method 'endswith' of 'str' objects}
        6    0.000    0.000    0.000    0.000 {method 'lower' of 'str' objects}
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/codecs.py:331(getstate)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:1495(parse_expression)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:353(fileno)
        4    0.000    0.000    0.000    0.000 {method 'fileno' of '_io.BufferedReader' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:95(_check_file_like)
        1    0.000    0.000    0.000    0.000 {method 'acquire' of '_thread.RLock' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/compat/py3k.py:49(isfileobj)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/tokenize.py:614(generate_tokens)
        2    0.000    0.000    0.000    0.000 {built-in method _struct.unpack}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/inspect.py:3177(bind)
        7    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/inspect.py:2695(kind)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:346(flush)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/fromnumeric.py:1862(_nonzero_dispatcher)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/sre_compile.py:477(_get_iscased)
        1    0.000    0.000    0.000    0.000 <__array_function__ internals>:177(atleast_1d)
        4    0.000    0.000    0.000    0.000 {method 'partition' of 'str' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/selection.py:237(__init__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/shape_base.py:218(_vhstack_dispatcher)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/memory.py:402(_format_hint)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/base.py:665(__init__)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/format.py:194(_check_version)
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_compression.py:160(tell)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/bz2.py:140(writable)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_compression.py:63(seekable)
        9    0.000    0.000    0.000    0.000 {method 'upper' of 'str' objects}
        1    0.000    0.000    0.000    0.000 {built-in method builtins.sum}
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_compression.py:36(readable)
        4    0.000    0.000    0.000    0.000 {built-in method builtins.iter}
        2    0.000    0.000    0.000    0.000 {method '__enter__' of '_io._IOBase' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:397(readline)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/_collections_abc.py:409(__subclasshook__)
        1    0.000    0.000    0.000    0.000 {method 'values' of 'mappingproxy' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/topology/base.py:117(__init__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:132(_resolve_formats)
        4    0.000    0.000    0.000    0.000 {built-in method _stat.S_ISLNK}
        5    0.000    0.000    0.000    0.000 {method 'reverse' of 'list' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:594(__init__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:526(__init__)
        2    0.000    0.000    0.000    0.000 {method 'popleft' of 'collections.deque' objects}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/lib/util.py:1060(asiterable)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/enum.py:801(value)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/inspect.py:2775(__init__)
        3    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/universe.py:479(universe)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/arraysetops.py:125(_unpack_tuple)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/chemfiles.py:171(_format_hint)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/multiarray.py:341(where)
        2    0.000    0.000    0.000    0.000 {method 'keys' of 'dict' objects}
        6    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/inspect.py:2683(name)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/logging/__init__.py:1307(disable)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:2556(n_atoms)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:2448(<genexpr>)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:368(readable)
        1    0.000    0.000    0.000    0.000 {method 'items' of 'mappingproxy' objects}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/inspect.py:3002(parameters)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/topologyattrs.py:559(__init__)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/core/groups.py:3227(<listcomp>)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/arraysetops.py:519(_in1d_dispatcher)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/multiarray.py:148(concatenate)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:371(writable)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/numeric.py:2190(_allclose_dispatcher)
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/contextlib.py:530(__enter__)
        1    0.000    0.000    0.000    0.000 {method 'release' of '_thread.RLock' objects}
        1    0.000    0.000    0.000    0.000 {method 'update' of 'dict' objects}
        2    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/gzip.py:374(seekable)
        2    0.000    0.000    0.000    0.000 {method 'pop' of 'collections.deque' objects}
        2    0.000    0.000    0.000    0.000 {method 'append' of 'collections.deque' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/twodim_base.py:438(_trilu_dispatcher)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/lib/arraysetops.py:133(_unique_dispatcher)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/numpy/core/shape_base.py:19(_atleast_1d_dispatcher)
        1    0.000    0.000    0.000    0.000 {method 'pop' of 'dict' objects}
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/base.py:601(close)
        1    0.000    0.000    0.000    0.000 /home/leonardo/.conda/envs/ge/lib/python3.10/site-packages/MDAnalysis/coordinates/base.py:605(__enter__)


