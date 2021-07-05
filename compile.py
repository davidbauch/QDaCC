# Extracts all runtimes from logfile

import os
import sys
from time import sleep
import subprocess
import multiprocessing
from functools import partial

def get_all_source_files(file_path="/", source_path="source", bin_path="obj", exclude_dirs=[], extension_source=".cpp", extension_obj=".o"):
    ret = []
    for root, dirs, files in os.walk(os.path.join(file_path, source_path)):
        for file in files:
            excluded = False
            for exdir in exclude_dirs:
                if exdir in root:
                    excluded = True
                    break
            if extension_source in file and not excluded:
                src_path = os.path.join(root, file).replace("\\", "/")
                obj_path = os.path.join(file_path, bin_path, file.replace(
                    extension_source, extension_obj)).replace("\\", "/")
                obj = file.replace(extension_source, extension_obj)
                ret.append([file, obj, src_path, obj_path])
    return ret


def compile_single(cset, force_recompile = False, compiler="g++", libs="", args="", cerr_to_file = False):
    name, obj, path_src, path_obj = cset[0], cset[1], cset[2], cset[3]
    cerr = "2>&1  | tee  error/{}_cerr.txt".format(name) if cerr_to_file else ""
    src_time = os.path.getmtime(path_src)
    obj_time = 0
    if not os.path.exists(path_obj):
        add_to_compile = True
    else:
        obj_time = os.path.getmtime(path_obj)
    try:
        if (obj_time < src_time) or force_recompile:
            print("Compiling {} ... ".format(name), end="\n", flush=True)
            cur = os.system('{} -c "{}" {} {} {}'.format(compiler, path_src, libs, args, cerr))
            #print("Moving object file ... ", end="", flush=True)
            os.replace(obj, path_obj)
            #print("Done!", flush=True)
            return "NEW"
        else:
            return "SKIPPED"
    except Exception as e:
        return "FAILED"
    return "FAILED"
        

def compile_all_object_files(data, bin_path="obj", extension_obj=".o", compiler="g++", libs="", args="", force_recompile=False, cerr_to_file=False, threads = -1):
    if not os.path.exists(bin_path):
        print("Bin path did not exist, creating '{}'".format(bin_path))
        os.makedirs(bin_path)
    if force_recompile:
        print("Removing old {} files ... ".format(extension_obj), end="")
        for file in os.listdir(bin_path):
            if file.endswith(extension_obj):
                os.remove(os.path.join(bin_path, file))
        print("Done!")
    if "-g" in args or "-g" in compiler:
        print("Compiling with debug information!")
    ret = []
    if cerr_to_file:
        os.makedirs("error", exist_ok=True)

    wrapper = partial(compile_single,cerr_to_file=cerr_to_file, compiler=compiler, libs=libs, args=args, force_recompile=force_recompile)
    with multiprocessing.Pool(processes=8) as pool:
        ret = pool.map(wrapper,[cset for cset in data])
    succesful = ret.count("FAILED") == 0
    if ret.count("FAILED") and ret.count("NEW") > 0:
        print("Compiling was {}.".format(
            "succesful." if succesful else "unsuccesful! ({} failed)".format(ret.count("FAILED"))))
    return succesful


def compile_main_program(path="", target="main.cpp", name="main.o", extension_obj=".o", bin_path="obj", compiler="g++", libs="", args="", copy_to="", add_base_path_to_obj=True):
    print("Compiling {} into {} ... ".format(target, name), end="", flush=True)
    final_bin_path = ""
    if not isinstance(bin_path, list):
        bin_path = list(bin_path)
    for binpath in bin_path:
        objpath = os.path.join(
            path, binpath) if add_base_path_to_obj else binpath
        objpath = objpath + "/*"+extension_obj
        objpath = objpath.replace("//", "/")
        final_bin_path += objpath + " "
    #os.system('{} "{}" -o "{}" "{}" {} {}'.format( compiler, os.path.join(path,target), os.path.join(path,name), final_bin_path, libs, args ))
    subprocess.Popen('{} "{}" -o "{}" {} {} {}'.format(compiler, os.path.join(path, target),
                                                       os.path.join(path, name), final_bin_path, libs, args), shell=True).wait()
    print("Done!", flush=True)
    if copy_to != "":
        print("Moving {} to {} ...".format(name, copy_to), end="")
        os.replace(name, copy_to)
        print("Done!")


if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    
    if "-frc" in sys.argv:
        if os.path.isfile(os.path.join(path,".fullcompile")):
            a = input("Last compile was a forced fullcompile. Press [ENTER] if you want to fullcompile again. Any other input will cancel the compilation.")
            if len(a) > 0:
                exit()
        else:
            with open(os.path.join(path,".fullcompile"),"w") as f:
                pass
    elif os.path.isfile(os.path.join(path,".fullcompile")):
        os.remove(os.path.join(path,".fullcompile"))
    
    platform = sys.platform
    if "-laptop" in sys.argv:
        platform = "win32"
    elif "-cluster" in sys.argv:
        platform = "cluster"
    print("Platform is {}".format(platform))

    msyspath = "msys2"

    # Get Program version
    version = "0.0"
    with open(os.path.join(path, "include/system/parameters.h")) as f:
        data = f.readlines()
        for line in data:
            if "GLOBAL_PROGRAM_VERSION" in line:
                version = line.split()[-1][1:-1]
    print("Will try to compile into version {}".format(version))

    force_recompile = True if "-frc" in sys.argv else False
    cerr_to_file = True if "-cerr" in sys.argv else False
    debug_info = " -g " if "-g" in sys.argv else ""

    bin_path = {'win32': "obj/" +
                ("LAPN" if "-laptop" in sys.argv else "win"), "darwin": "obj/MAC", "cluster" : "obj/CLUSTER"}
    compiler = {'win32': "g++", 'darwin': "g++-9", "cluster" : "g++"}

    f = get_all_source_files(
        path, exclude_dirs=["ALGLIB"], bin_path=bin_path[platform])

    # Debug : -g -O0, else: -O3

    copy_to = {'win32': "../../Threadhandler/QDLC-{}.exe".format(
        version) if "-th" in sys.argv else "", 'darwin': "/Users/davidbauch/bin/QDLC-{}.out".format(version), "cluster" : "/"}
    libs_obj = {'win32': "-Wno-volatile -std=c++2a -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs",
                'darwin': "-Wno-volatile -std=c++17 -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs",
                "cluster" : "-Wno-volatile -std=c++17 -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs"}
    include_obj = {'win32': '-I"C:\{}\myinclude" -I"C:/Users/david/OneDrive - Universität Paderborn/Kot/BP/QDLC-C/include" -I"C:/Users/david/OneDrive - Universität Paderborn/Kot/BP/QDLC-C/external"'.format(msyspath), 
                    'darwin': "-I'/Users/davidbauch/OneDrive - Universität Paderborn/Kot/BP/QDLC-C/include' -I'/Users/davidbauch/OneDrive - Universität Paderborn/Kot/BP/QDLC-C/external'",
                    "cluster" : "-I'/include' -I'/external'"}

    libs_final = {'win32': '-Wno-volatile -std=c++2a -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs {}'.format("-static" if "-static" in sys.argv else ""), 
                    'darwin': "-Wno-volatile -std=c++17 -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs {}".format("-static" if "-static" in sys.argv else ""),
                    "cluster" : "-Wno-volatile -std=c++17 -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs {}".format("-static" if "-static" in sys.argv else "")}
    bin_final = {'win32': ["obj/"+("LAPN" if "-laptop" in sys.argv else "win"), "external/ALGLIB/"+("LAPN" if "-laptop" in sys.argv else "WIN")], 
                'darwin': ["obj/MAC", "external/ALGLIB/MAC"],
                "cluster" : ["obj/CLUSTER", "external/ALGLIB/CLUSTER"]}
    succesful = compile_all_object_files(f, libs=libs_obj[platform], args=include_obj[platform], force_recompile=force_recompile,
                                         bin_path=bin_path[platform], compiler=compiler[platform]+debug_info, cerr_to_file=cerr_to_file)
    # if succesful:
    add_basepath = {'win32': False, 'darwin': False, "cluster" : False}
    compile_main_program(path=path, libs=libs_final[platform], args=include_obj[platform], copy_to=copy_to[platform],
                         bin_path=bin_final[platform], compiler=compiler[platform]+debug_info, add_base_path_to_obj=add_basepath[platform])
