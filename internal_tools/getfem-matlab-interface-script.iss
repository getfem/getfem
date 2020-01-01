#define MyAppName "Getfem-5.2-Matlab2016"
#define MyAppVersion "5.2"
#define MyAppPublisher "Getfem project"
#define MyAppURL "http://www.getfem.org"

[Setup]
AppId={{ACCE7B33-300B-4CAB-9D51-EAA9E007A1C0}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
;AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}
DefaultDirName={pf}\{#MyAppName}
DefaultGroupName={#MyAppName}
LicenseFile=C:\msys\home\yves\getfem-5.1\COPYING
OutputDir=C:\msys\home\yves
OutputBaseFilename=getfem5.2-matlab2016-interface-setup
Compression=lzma
SolidCompression=yes
ChangesEnvironment=yes
DisableDirPage=no

[Code]
procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
var
  Path, AppDir: string;
  Index: Integer;
begin
  if CurUninstallStep = usUninstall then
  begin
    if RegQueryStringValue(HKEY_LOCAL_MACHINE,
      'SYSTEM\CurrentControlSet\Control\Session Manager\Environment\',
      'MATLABPATH', Path) then
    begin
      AppDir := ExpandConstant('{app}');
      Index := Pos(AppDir, Path);
      Delete(Path, Index-1, Length(AppDir)+1);
      RegWriteStringValue(HKEY_LOCAL_MACHINE,
        'SYSTEM\CurrentControlSet\Control\Session Manager\Environment\',
        'MATLABPATH', Path);
    end;
  end;
end;

[Code]
function NeedsAddPath(Param: string): boolean;
var
  OrigPath: string;
begin
  if not RegQueryStringValue(HKEY_LOCAL_MACHINE,
    'SYSTEM\CurrentControlSet\Control\Session Manager\Environment',
    'MATLABPATH', OrigPath)
  then begin
    Result := True;
    exit;
  end;
  { look for the path with leading and trailing semicolon }
  { Pos() returns 0 if not found }
  Result := Pos(';' + Param + ';', ';' + OrigPath + ';') = 0;
end;

[Registry]
; Disabled because set MATLABPATH disturb Matlab 2016 ...
; Root: HKLM; SubKey: "SYSTEM\CurrentControlSet\Control\Session Manager\Environment\"; ValueType: string; ValueName: "MATLABPATH"; ValueData: "{reg:HKLM\SYSTEM\CurrentControlSet\Control\Session Manager\Environment\,MATLABPATH};{app}"; Check: NeedsAddPath('{app}'); Flags: uninsdeletekeyifempty

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "armenian"; MessagesFile: "compiler:Languages\Armenian.islu"
Name: "brazilianportuguese"; MessagesFile: "compiler:Languages\BrazilianPortuguese.isl"
Name: "catalan"; MessagesFile: "compiler:Languages\Catalan.isl"
Name: "corsican"; MessagesFile: "compiler:Languages\Corsican.isl"
Name: "czech"; MessagesFile: "compiler:Languages\Czech.isl"
Name: "danish"; MessagesFile: "compiler:Languages\Danish.isl"
Name: "dutch"; MessagesFile: "compiler:Languages\Dutch.isl"
Name: "finnish"; MessagesFile: "compiler:Languages\Finnish.isl"
Name: "french"; MessagesFile: "compiler:Languages\French.isl"
Name: "german"; MessagesFile: "compiler:Languages\German.isl"
Name: "greek"; MessagesFile: "compiler:Languages\Greek.isl"
Name: "hebrew"; MessagesFile: "compiler:Languages\Hebrew.isl"
Name: "hungarian"; MessagesFile: "compiler:Languages\Hungarian.isl"
Name: "italian"; MessagesFile: "compiler:Languages\Italian.isl"
Name: "japanese"; MessagesFile: "compiler:Languages\Japanese.isl"
Name: "nepali"; MessagesFile: "compiler:Languages\Nepali.islu"
Name: "norwegian"; MessagesFile: "compiler:Languages\Norwegian.isl"
Name: "polish"; MessagesFile: "compiler:Languages\Polish.isl"
Name: "portuguese"; MessagesFile: "compiler:Languages\Portuguese.isl"
Name: "russian"; MessagesFile: "compiler:Languages\Russian.isl"
Name: "scottishgaelic"; MessagesFile: "compiler:Languages\ScottishGaelic.isl"
Name: "serbiancyrillic"; MessagesFile: "compiler:Languages\SerbianCyrillic.isl"
Name: "serbianlatin"; MessagesFile: "compiler:Languages\SerbianLatin.isl"
Name: "slovenian"; MessagesFile: "compiler:Languages\Slovenian.isl"
Name: "spanish"; MessagesFile: "compiler:Languages\Spanish.isl"
Name: "turkish"; MessagesFile: "compiler:Languages\Turkish.isl"
Name: "ukrainian"; MessagesFile: "compiler:Languages\Ukrainian.isl"

[Files]
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\*.m";                   DestDir: "{app}";                   Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\*.mex*";                DestDir: "{app}";                   Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\private\*.m";           DestDir: "{app}\private";           Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfContStruct\*.m";     DestDir: "{app}\@gfContStruct";     Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfCvStruct\*.m";       DestDir: "{app}\@gfCvStruct";       Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfEltm\*.m";           DestDir: "{app}\@gfEltm";           Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfFem\*.m";            DestDir: "{app}\@gfFem";            Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfGeoTrans\*.m";       DestDir: "{app}\@gfGeoTrans";       Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfGlobalFunction\*.m"; DestDir: "{app}\@gfGlobalFunction"; Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfInteg\*.m";          DestDir: "{app}\@gfInteg";          Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfLevelSet\*.m";       DestDir: "{app}\@gfLevelSet";       Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfMesh\*.m";           DestDir: "{app}\@gfMesh";           Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfMeshFem\*.m";        DestDir: "{app}\@gfMeshFem";        Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfMeshIm\*.m";         DestDir: "{app}\@gfMeshIm";         Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfMeshImData\*.m";     DestDir: "{app}\@gfMeshImData";     Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfMeshLevelSet\*.m";   DestDir: "{app}\@gfMeshLevelSet";   Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfMesherObject\*.m";   DestDir: "{app}\@gfMesherObject";   Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfModel\*.m";          DestDir: "{app}\@gfModel";          Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfPrecond\*.m";        DestDir: "{app}\@gfPrecond";        Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfSlice\*.m";          DestDir: "{app}\@gfSlice";          Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\src\matlab\@gfSpmat\*.m";          DestDir: "{app}\@gfSpmat";          Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\tests\matlab\*.m";                 DestDir: "{app}\tests";             Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\tests\matlab\private\*.m";         DestDir: "{app}\tests\private";     Flags: ignoreversion
Source: "C:\msys\home\yves\getfem-5.2\interface\tests\meshes\*.msh";               DestDir: "{app}\meshes";            Flags: ignoreversion


