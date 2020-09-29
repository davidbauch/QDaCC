# Extracts all runtimes from logfile

import os
import sys
from time import sleep

def get_all_source_files( file_path = "/", source_path = "source", bin_path = "obj", exclude_dirs = [], extension_source = ".cpp", extension_obj = ".o" ):
    ret = []
    for root, dirs, files in os.walk(os.path.join(file_path, source_path)):
        for file in files:
            excluded = False
            for exdir in exclude_dirs:
                if exdir in root:
                    excluded = True
                    break
            if extension_source in file and not excluded:
                src_path = os.path.join(root,file).replace("\\","/")
                obj_path = os.path.join(file_path,bin_path,file.replace(extension_source,extension_obj)).replace("\\","/")
                obj = file.replace(extension_source,extension_obj)
                ret.append([file, obj, src_path, obj_path])
    return ret

def compile_all_object_files( data, bin_path = "obj", extension_obj = ".o", compiler = "g++", libs = "", args = "", force_recompile = False ):
    if force_recompile:
        print("Removing old {} files ... ".format(extension_obj),end="")
        for file in os.listdir(bin_path):
            if file.endswith(extension_obj):
                os.remove(os.path.join(bin_path,file))
        print("Done!")
    ret = []
    succesful = True
    for set in data:
        name,obj,path_src,path_obj = set[0], set[1], set[2], set[3]
        add_to_compile = False
        src_time = os.path.getmtime( path_src )
        obj_time = 0
        if not os.path.exists( path_obj ):
            add_to_compile = True
        else:
            obj_time = os.path.getmtime( path_obj )
        if ( obj_time < src_time ) or force_recompile:
            print("Compiling {} ... ".format(name),end="",flush=True)
            cur = os.system('{} -c "{}" {} {}'.format(compiler,path_src,libs,args))
            succesful = succesful and cur
            print("Moving object file ... ",end="",flush=True)
            os.replace(obj,path_obj)
            print("Done!",flush=True)
    if len(ret) > 0:
        print("Compiling was {}.".format("succesful." if succesful else "unsuccesful!"))
    return succesful

def compile_main_program( path = "", target = "main.cpp", name = "main.exe", extension_obj = ".o", bin_path = "obj", compiler = "g++", libs = "", args = "", copy_to = ""):
    print("Compiling {} into {} ... ".format(target,name),end="",flush=True)
    os.system('{} "{}" -o "{}" "{}/*{}" {} {}'.format( compiler, os.path.join(path,target), os.path.join(path,name), os.path.join(path,bin_path), extension_obj, libs, args ))
    print("Done!",flush=True)
    if copy_to != "":
        print("Moving {} to {} ...".format(name,copy_to),end="")
        os.replace(name,copy_to)
        print("Done!")



if __name__ == "__main__":
    path = os.path.dirname(os.path.realpath(__file__))
    copy_to = "../../Threadhandler/QDLC-2.0.exe" if "-th" in sys.argv else ""
    force_recompile = True if "-frc" in sys.argv else False
    f = get_all_source_files(path, exclude_dirs=["ALGLIB"])
    succesful = compile_all_object_files(f, libs = "-std=c++2a -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs", args = '-I"C:\msys2\myinclude" -I"C:/Users/david/OneDrive - Universität Paderborn/Kot/BP/QDLC-C/include"', force_recompile=force_recompile)
    #if succesful:
    compile_main_program(path=path, libs = '"C:/Users/david/OneDrive - Universität Paderborn/Kot/BP/QDLC-C/obj/ALGLIB/WIN/*.o" -std=c++2a -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs', args = '-I"C:\msys2\myinclude" -I"C:/Users/david/OneDrive - Universität Paderborn/Kot/BP/QDLC-C/include"', copy_to=copy_to)
