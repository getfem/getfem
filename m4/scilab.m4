AC_DEFUN([AC_CHECK_SCILAB],
[
  AH_TEMPLATE(HAVE_SCILAB,
   [Defined to 1 if Scilab is present on the system])

  AC_ARG_ENABLE(scilab,
  [AS_HELP_STRING([--enable-scilab],[turn on/off scilab support])],
  [case "${enableval}" in
  	 yes) usescilab=YES ;;
   	 no)  usescilab=NO ;;
   	 *) AC_MSG_ERROR([bad value ${enableval} for --enable-scilab]) ;;
   esac],[usescilab=NO])

  AC_ARG_WITH(scilab_prefix,
		AC_HELP_STRING([--with-scilab-prefix=DIR],[Set the path to Scilab]),
		[with_scilab_prefix=$withval],
		[with_scilab_prefix='yes']
		)

  AC_ARG_WITH(scilab_version,
		AC_HELP_STRING([--with-scilab-version="major.minor.micro"],[Set the required Scilab version]),
		[with_scilab_version=$withval],
		[with_scilab_version='yes']
		)

  AC_ARG_WITH(scilab_toolbox_dir,
		AC_HELP_STRING([--with-scilab-toolbox-dir=DIR],[Set the path to the toolbox installation directory]),
		[with_scilab_toolbox_dir=$withval],
		[with_scilab_toolbox_dir='yes']
		)

  if test "x$usescilab" == "xYES"
  then
    dnl if version is not set by the user, we get the current scilab version as required version
    if test -z $REQUIRED_SCILAB_MAJOR
    then
      REQUIRED_SCILAB_MAJOR=`echo "$SCILAB_VERSION" | sed  "s/.*\([[0-9]]\+\)[[.]]\([[0-9]]\+\)[[.]]\([[0-9]]\+\)/\1/"`
    fi
    if test -z $REQUIRED_SCILAB_MINOR
    then
      REQUIRED_SCILAB_MINOR=`echo "$SCILAB_VERSION" | sed  "s/.*\([[0-9]]\+\)[[.]]\([[0-9]]\+\)[[.]]\([[0-9]]\+\)/\2/"`
    fi
    if test -z $REQUIRED_SCILAB_MICRO
    then
      REQUIRED_SCILAB_MICRO=`echo "$SCILAB_VERSION" | sed  "s/.*\([[0-9]]\+\)[[.]]\([[0-9]]\+\)[[.]]\([[0-9]]\+\)/\3/"`
    fi

    dnl check for Scilab

    if test "x$with_scilab_prefix" != "xyes"
    then
      if test -x "$with_scilab_prefix/bin/scilab"
      then 
        AC_MSG_RESULT([Scilab binary program was found in $with_scilab_prefix])
      else
        AC_MSG_ERROR([Scilab binary program was not found in $with_scilab_prefix/bin])
      fi
      SCILAB_EXE="$with_scilab_prefix/bin/scilab"
      AC_DEFINE(HAVE_SCILAB)
    else
      AC_CHECK_PROG([has_scilab],[scilab],[yes],[no], , )
      if test x$has_scilab = xno; then
        AC_MSG_ERROR([[Scilab binary program was found in your PATH], your PATH is $PATH])
      fi
      SCILAB_EXE="scilab"
      AC_DEFINE(HAVE_SCILAB)
    fi

    dnl test if SCI is defined. If not, test if scilab is in path,
    dnl then get SCI from scilab -nwni.
    if test -z "$SCI"; then
      cmd='F=mopen("getpath.incl","w");
           mfprintf(F,SCI);
           mclose(F);exit;'
      echo "$cmd" > getpath.sci
      $SCILAB_EXE -nw -f getpath.sci >/dev/null
      SCILAB_DIR=`cat getpath.incl`
      rm -f getpath.sci getpath.incl 
    else
      SCILAB_DIR="$SCI"
    fi

    dnl get scilab version
    cmd='F=mopen("version.incl","w");
         ver=getversion();
         mfprintf(F,ver);
         mclose(F);exit;'
    echo "$cmd" > version.sci
    $SCILAB_EXE -nwni -f version.sci >/dev/null
    SCILAB_VERSION=`cat version.incl`
    rm -f version.sci version.incl 

    scilab_tmp_version=`echo $SCILAB_VERSION | sed -r "s/.*(branch).*/\1/"`

    if test "x$scilab_tmp_version" = "xbranch"
    then
      SCILAB_VERSION_MAJOR=-1
      SCILAB_VERSION_MINOR=-1
      SCILAB_VERSION_MICRO=-1
    else
      SCILAB_VERSION_MAJOR=`echo "$SCILAB_VERSION" | sed -r "s/.*([[0-9]]+)[[.]]([[0-9]]+)[[.]]([[0-9]]+)/\1/"`
      SCILAB_VERSION_MINOR=`echo "$SCILAB_VERSION" | sed -r "s/.*([[0-9]]+)[[.]]([[0-9]]+)[[.]]([[0-9]]+)/\2/"`
      SCILAB_VERSION_MICRO=`echo "$SCILAB_VERSION" | sed -r "s/.*([[0-9]]+)[[.]]([[0-9]]+)[[.]]([[0-9]]+)/\3/"`

      if test $SCILAB_VERSION_MAJOR -lt $REQUIRED_SCILAB_MAJOR
      then
        AC_MSG_ERROR([scilab major version does not match])
      else
        if test $SCILAB_VERSION_MINOR -lt $REQUIRED_SCILAB_MINOR
        then
          AC_MSG_ERROR([scilab minor version does not match])
        else
          if test $SCILAB_VERSION_MICRO -lt $REQUIRED_SCILAB_MICRO
          then
            AC_MSG_ERROR([scilab micro version does not match])
          fi
        fi
      fi
    fi

    if test "xwith_scilab_toolbox_dir" != "xyes"
    then
      SCILAB_TOOLBOX_DIR="$with_scilab_toolbox_dir"
    else
      SCILAB_TOOLBOX_DIR="$SCILAB_DIR/contrib/$PACKAGE_NAME-$PACKAGE_VERSION"
    fi
  fi

  AM_CONDITIONAL(BUILDSCILAB, test x$usescilab = xYES)

  AC_SUBST(HAVE_SCILAB)
  AC_SUBST(SCILAB_EXE)
  AC_SUBST(SCILAB_DIR)
  AC_SUBST(SCILAB_TOOLBOX_DIR)
  AC_SUBST(SCILAB_VERSION_MAJOR)
  AC_SUBST(SCILAB_VERSION_MINOR)
  AC_SUBST(SCILAB_VERSION_MICRO)
  AC_SUBST(SCILAB_VERSION)
])
