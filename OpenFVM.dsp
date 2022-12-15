# Microsoft Developer Studio Project File - Name="OpenFVM" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=OpenFVM - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "OpenFVM.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "OpenFVM.mak" CFG="OpenFVM - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "OpenFVM - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "OpenFVM - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "OpenFVM - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "examples"
# PROP Intermediate_Dir "debug"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x1009 /d "NDEBUG"
# ADD RSC /l 0x1009 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "OpenFVM - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "examples"
# PROP Intermediate_Dir "serial/debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /Gm /GX /ZI /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x1009 /d "_DEBUG"
# ADD RSC /l 0x1009 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "OpenFVM - Win32 Release"
# Name "OpenFVM - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\serial\source\bcond.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\decomp.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\fill.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\gamma.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\geocalc.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\gradient.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\ioutils.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\main.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\material.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\mesh.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\msolver.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\octree.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\param.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\parser.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\post.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\pressure.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\rcm.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\reorder.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\restart.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\setup.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\solve.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\temperature.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\ttime.c
# End Source File
# Begin Source File

SOURCE=.\serial\source\velocity.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\serial\source\bcond.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\decomp.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\fill.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\gamma.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\geocalc.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\globals.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\gradient.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\ioutils.h
# End Source File
# Begin Source File

SOURCE=.\serial\laspack\laspack.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\material.h
# End Source File
# Begin Source File

SOURCE=.\serial\laspack\matrix.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\mesh.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\msolver.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\octree.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\param.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\parser.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\post.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\pressure.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\rcm.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\reorder.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\restart.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\setup.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\solve.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\stress.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\temperature.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\ttime.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\variables.h
# End Source File
# Begin Source File

SOURCE=.\serial\source\velocity.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\serial\laspack\laspack.lib
# End Source File
# End Target
# End Project
