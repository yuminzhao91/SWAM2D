\NeedsTeXFormat{LaTeX2e}[1996/12/01]
\ProvidesClass{gnulike}[2017/01/10 v1.0]
%%
%% Pass all options directly to 'article'
%%
%\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
%\ProcessOptions\relax
%%
%% prefered fontsize, papersize, and type
%%
\LoadClass[oneside,11pt,a4paper]{report}
%%
%% prefered fonts: Times New Roman becomes the standard font (text)
%%
\RequirePackage[T1]{fontenc}
\RequirePackage{times}
\RequirePackage{mathptmx}
\RequirePackage{titlesec}
%\RequirePackage[hidelinks,colorlinks=true,citecolor=vibrisblue2,linkcolor=vibrisblue2,urlcolor=blue]{hyperref}
\RequirePackage[hidelinks]{hyperref}
%%
%% Page size formats
%%
\RequirePackage[a4paper]{geometry}
\geometry{left=2.5cm,right=2.5cm,top=2.5cm,bottom=2.5cm}

\RequirePackage{color}

% >> Tango color palette
% <http://tango.freedesktop.org/static/cvs/tango-art-tools/palettes/Tango-Palette.png>

% yellow
\definecolor{butter1}{RGB}{252,233,79}
\definecolor{butter2}{RGB}{237,212,0}
\definecolor{butter3}{RGB}{196,160,0}

% orange
\definecolor{orange1}{RGB}{252,175,62}
\definecolor{orange2}{RGB}{245,121,0}
\definecolor{orange3}{RGB}{206,92,0}
\definecolor{orange1t}{RGB}{252,200,87}

% brown
\definecolor{chocolate1}{RGB}{233,185,110}
\definecolor{chocolate2}{RGB}{193,125,17}
\definecolor{chocolate3}{RGB}{143,89,2}

% green
\definecolor{chameleon1}{RGB}{138,226,52}
\definecolor{chameleon2}{RGB}{115,210,22}
\definecolor{chameleon3}{RGB}{78,154,6}

% blue
\definecolor{skyblue1}{RGB}{114,159,207}
\definecolor{skyblue2}{RGB}{52,101,164}
\definecolor{skyblue3}{RGB}{32,74,135}

% purple
\definecolor{plum1}{RGB}{173,127,168}
\definecolor{plum2}{RGB}{117,80,123}
\definecolor{plum3}{RGB}{92,53,102}

% red
\definecolor{scarlet1}{RGB}{239,41,41}
\definecolor{scarlet2}{RGB}{204,0,0}
\definecolor{scarlet3}{RGB}{164,0,0}

% gray
\definecolor{aluminium1}{RGB}{238,238,236}
\definecolor{aluminium2}{RGB}{211,215,207}
\definecolor{aluminium3}{RGB}{186,189,182}
\definecolor{aluminium4}{RGB}{136,138,133}
\definecolor{aluminium5}{RGB}{85,87,83}
\definecolor{aluminium6}{RGB}{46,52,54}

%%
%% bibliography
%%
\RequirePackage{natbib}
\gdef\bibfont{\fontsize{10}{10}\selectfont}
\setlength{\bibhang}{0em}
\bibliographystyle{apalike}

%%
%% set header and footer
%%
\RequirePackage{fancyhdr}
\RequirePackage{graphicx}
\pagestyle{fancy}
\fancyhf{}
\rhead{\textcolor{aluminium4}{\thepage}}
\lhead{\textcolor{aluminium4}{\leftmark}}
\renewcommand{\headrulewidth}{1pt}% 2pt header rule
\renewcommand{\headrule}{\hbox to\headwidth{%
  \color{aluminium4}\leaders\hrule height \headrulewidth\hfill}}

\fancypagestyle{plain}{%
  \fancyhf{}%
  \rhead{\textcolor{aluminium4}{\thepage}}
  \lhead{\textcolor{aluminium4}{\leftmark}}
}

\fancypagestyle{contents}{%
  \fancyhf{}%
  \rhead{}
  \lhead{}
}
