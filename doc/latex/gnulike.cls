%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is file 'stp2012.cls' %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Made by: Bobby Hak
%%
%% This class has been tested with:
%%     MikTeX version 2.9
%%
%% It is based on the standard LaTeX2e "article" class and some
%% standard packages:
%%    times    - for the font Times New Roman, 11pt
%%    mathptmx - for the font Times New Romen with math support
%%    fancyhdr - for defining headnote and footnote
%%    graphicx - for including logo
%%    geometry - for putting page dimensions correct
%%    natbib   - for referencing
%%
%% All packages should be availble in any normal distribution of LaTeX,
%% or can be downloaded from the CTAN homepage:
%%    www.ctan.org
%%
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
%\geometry{
% a4paper,
% left=20mm,
% right=20mm,
% }
\geometry{left=2.5cm,right=2.5cm,top=2.5cm,bottom=2.5cm}
%\geometry{head=0.68in,headsep=0.09in,foot=0.5in}
%\addtolength{\footskip}{-1.0\baselineskip}
%%
%% some others
%%
%\setcounter{secnumdepth}{0} % no numbering of sections
%\setlength{\parindent}{0pt} % no indent for new paragraphs
%\setlength{\parskip}{11pt} % skip one line between paragraphs
%\setlength{\partopsep}{0pt} % skip between environments and paragraphs
%\setlength{\topsep}{0pt} % skip between environments and paragraphs
%%
%% bibliography
%%
\RequirePackage{natbib} % special citing commands
\gdef\bibfont{\fontsize{10}{10}\selectfont}
%\newcommand{\bibfont}{\fontsize{10}{10}\selectfont} % bibliography fontsize
%\setlength{\bibsep}{1pt} % no lines between bibitems
\setlength{\bibhang}{0em} % indentation of second line in bibitem
\bibliographystyle{apalike}  %{firstbreak} % stylefile
%%
%% set header and footer
%%
\RequirePackage{fancyhdr}
\RequirePackage{graphicx}
\pagestyle{fancy}
\fancyhf{}
\rhead{\thepage}
\lhead{\leftmark}
%\rfoot{Page \thepage}
%\pagestyle{fancy} % load package
%\fancyhf{} % clear styles
%\fancyhead{test} % input logo top right, even and odd side
%\fancyfoot[C]{\vspace{-1pt}\fontsize{10}{10}\selectfont{%
%\color{vibrisblue2}VIBRIS Project -- Intermediate Report\\
%Physical and Numerical Experiments -- Damien Pageot
%}} % footnote
\renewcommand{\headrulewidth}{0pt} % no line in header

\fancypagestyle{plain}{%
  \fancyhf{}%
  \rhead{\thepage}
  \lhead{\nouppercase{\leftmark}}
  %\fancyfoot[C]{\thepage\ of \pageref{LastPage}}%
  %\renewcommand{\headrulewidth}{0pt}% Line at the header invisible
  %\renewcommand{\footrulewidth}{0.4pt}% Line at the footer visible
}

\fancypagestyle{contents}{%
  \fancyhf{}%
  \rhead{}
  \lhead{}
  %\fancyfoot[C]{\thepage\ of \pageref{LastPage}}%
  %\renewcommand{\headrulewidth}{0pt}% Line at the header invisible
  %\renewcommand{\footrulewidth}{0.4pt}% Line at the footer visible
}

%\renewcommand{\footrulewidth}{0.4pt} % add line in footnote
%\pagenumbering{gobble} % define page numbering
%%
%% section headers and paragraphs
%%
%\renewcommand{\section}{\@startsection{section}%
%{1}% level
%{0pt}% indent
%{.1ex plus -.1ex}% beforeskip
%{.1ex plus -.1ex}% afterskip
%{\fontsize{11}{11}\selectfont\bf}} % font
%
%\renewcommand{\subsection}{\@startsection{subsection}%
%{2}% level
%{0pt}% indent
%{.1ex plus -.1ex}% beforeskip
%{.1ex plus -.1ex}% afterskip
%{\fontsize{11}{11}\selectfont\emph\bf}} % font
%%
%% caption font
%%
%\long\def\@makecaption#1#2{%
%  \vskip\abovecaptionskip
%  \sbox\@tempboxa{{\bf\itshape #1} \itshape #2}%
%  \ifdim \wd\@tempboxa >\hsize
%    {{\bf\itshape #1} \itshape #2\par}
%  \else
%    \hbox to\hsize{\hfil\box\@tempboxa\hfil}%
%  \fi
%  \vskip\belowcaptionskip}
% EOF

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
