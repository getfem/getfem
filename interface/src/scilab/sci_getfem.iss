;##############################################################################################################
; Inno Setup Install script for SciGetFem Module
; http://www.jrsoftware.org/isinfo.php
; Yann COLLETTE
; This file is released into the public domain
;##############################################################################################################
; modify this path where is sciglpk directory
#define BinariesSourcePath "E:\Scilab\toolboxes\getfem\interface\src\scilab"
;
#define SCIGETFEM_Module_version "3914"
#define CurrentYear "2011"
#define SCIGETFEM_ModuleDirFilename "scigetfem-rev3914"
;##############################################################################################################
[Setup]
; Debut Données de base à renseigner suivant version
SourceDir={#BinariesSourcePath}
AppName=SciGetFem Module
AppVerName=SciGetFem Module version rev3914
DefaultDirName={pf}/{#SCIGETFEM_ModuleDirFilename}
InfoAfterfile=readme.txt
LicenseFile=license.txt
WindowVisible=true
AppPublisher=Home
BackColorDirection=lefttoright
AppCopyright=Copyright © {#CurrentYear}
Compression=lzma/max
InternalCompressLevel=normal
SolidCompression=true
VersionInfoVersion={#SCIGETFEM_Module_version}
VersionInfoCompany=Home
;##############################################################################################################
[Files]
; Add here files that you want to add
Source: loader.sce; DestDir: {app}
Source: license.txt; DestDir: {app}
Source: readme.txt; DestDir: {app}
Source: etc\sci_getfem.quit; DestDir: {app}\etc
Source: etc\sci_getfem.start; DestDir: {app}\etc
Source: sci_gateway\loader_gateway.sce; DestDir: {app}\sci_gateway
Source: sci_gateway\c\loader.sce; DestDir: {app}\sci_gateway\c
Source: sci_gateway\c\*.dll; DestDir: {app}\sci_gateway\c
Source: macros\*.*; DestDir: {app}\macros;
Source: macros\overload\*.*; DestDir: {app}\macros\overload;
Source: demos\*.*; DestDir: {app}\demos;
Source: demos\data\*.*; DestDir: {app}\demos\data;
Source: jar\*.*; DestDir: {app}\jar
Source: src\c\loader.sce; DestDir: {app}\src\c
Source: src\c\*.dll; DestDir: {app}\src\c
;
;##############################################################################################################
