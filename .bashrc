# $Id$
#
# This file is sourced for shells which are interactive but not a
# login shell; however, since it is also sourced within
# $HOME/.bash_profile, the end result is that this file is sourced for all
# shells.
PATH=$PATH:/u/dorislee/Programs/Montage_v3.3/bin:/u/dorislee/Programs/athena4.2/bin



# Source global definitions
if [ -f /etc/bashrc ]; then
  . /etc/bashrc
fi

# User specific aliases and functions
if [ -f $HOME/.bash_aliases ]; then
  . $HOME/.bash_aliases
fi

#
# Environment variables
#
[[ "$PS1" ]] && {
  if [ -z "$DARWIN" ] ; then
    export PS1='${debian_chroot:+($debian_chroot)}\[\e[01;`case $EUID in 0) echo -n 31; ;; *) echo -n 32; ;; esac;`m\]\h\[\e[00m\]:\[\e[01;34m\]\w\[\e[00m\]\$ '
  else
    export PS1='${debian_chroot:+($debian_chroot)}\[\e[01;`case $EUID in 0) echo "31\c"; ;; *) echo -n 32; ;; esac;`m\]\h\[\e[00m\]:\[\e[01;34m\]\w\[\e[00m\]\$ '
  fi
}
export PAGER=less
export EDITOR=vim

# If this is an xterm set the title to user@host:dir
case $TERM in
  xterm*|rxvt*)
    PROMPT_COMMAND='echo -ne "\033]0;${USER}@${HOSTNAME}: ${PWD/$HOME/~}\007"'
    ;;
  *)
    ;;
esac

#
# Paths
#
# Standard-ish paths; some of these may be in place already, but if
# they're not they should probably go near the front of the pack.
for P in /sbin /usr/sbin /usr/local/bin /usr/local/sbin \
  /opt/local/bin /opt/local/sbin /usr/X11R6/bin; do
  if [ -d $P ]; then
    [[ $PATH =~ (^|:)$P($|:) ]] || \
    export PATH=${P}${PATH:+:$PATH}
  fi
done

for M in /usr/local/share/man /usr/local/man /opt/local/share/man ; do
  if [ -d $M ]; then
    [[ $MANPATH =~ (^|:)$M($|:) ]] || \
    export MANPATH=${M}${MANPATH:+:$MANPATH}
  fi
done

# Local Perl install
if [ -d $HOME/perl ]; then
  if [ -d $HOME/perl/lib/perl5 ]; then
    [[ $PERL5LIB =~ (^|:)$HOME/perl/lib/perl5($|:) ]] || \
    export PERL5LIB=$HOME/perl/lib/perl5${PERL5LIB:+:$PERL5LIB}
  fi

  if [ -d $HOME/perl/lib64/perl5 ]; then
    [[ $PERL5LIB =~ (^|:)$HOME/perl/lib64/perl5($|:) ]] || \
    export PERL5LIB=$HOME/perl/lib64/perl5${PERL5LIB:+:$PERL5LIB}
  fi

  if [ -d $HOME/perl/lib/perl5/site_perl ]; then
    [[ $PERL5LIB =~ (^|:)$HOME/perl/lib/perl5/site_perl($|:) ]] || \
    export PERL5LIB=$HOME/perl/lib/perl5/site_perl${PERL5LIB:+:$PERL5LIB}
  fi

  if [ -d $HOME/perl/lib64/perl5/site_perl ]; then
    [[ $PERL5LIB =~ (^|:)$HOME/perl/lib64/perl5/site_perl($|:) ]] || \
    export PERL5LIB=$HOME/perl/lib64/perl5/site_perl${PERL5LIB:+:$PERL5LIB}
  fi

  if [ -d $HOME/perl/share/man ]; then
    [[ $MANPATH =~ (^|:)$HOME/perl/share/man($|:) ]] || \
    export MANPATH=$HOME/perl/share/man:${MANPATH:+$MANPATH}
  fi

  if [ -d $HOME/perl/bin ]; then
    [[ $PATH =~ (^|:)$HOME/perl/bin($|:) ]] || \
    export PATH=$HOME/perl/bin${PATH:+:$PATH}
  fi
fi

# Local install paths - install things to $HOME/Installs directories, and
# they will automatically get the proper paths added.  Most things
# that use autoconf will do this with
# './configure --prefix=$HOME/Installs'
for D in $HOME/Installs/* ; do
  if [ -d $D/bin ]; then 
    [[ $PATH =~ (^|:)$D/bin($|:) ]] || \
    export PATH=$D/bin${PATH:+:$PATH}
  fi

  if [ -d $D/share/man ]; then 
    [[ $MANPATH =~ (^|:)$D/share/man($|:) ]] || \
    export MANPATH=$D/share/man:${MANPATH:+$MANPATH}
  fi
  
  if [ -d $D/man ]; then 
    [[ $MANPATH =~ (^|:)$D/man($|:) ]] || \
    export MANPATH=$D/man:${MANPATH:+$MANPATH}
  fi
done

# Add-ons for Condor
if [ -f /u/condor/condor-setup.sh ]; then
  [[ $PATH =~ (^|:)/u/condor/hosts/`hostname -s`/sbin($|:) ]] || {
    export PATH=${PATH:+$PATH:}/u/condor/hosts/`hostname -s`/sbin
    . /u/condor/condor-setup.sh
  }
fi

# Now make sure $HOME/bin is top of the list (or at least present)
if [ -d $HOME/bin ]; then
  [[ $PATH =~ (^|:)$HOME/bin($|:) ]] || \
  export PATH=$HOME/bin${PATH:+:$PATH}
fi

# Last, make sure to end MANPATH with a ':' to force include of system paths
# (some versions of man take care of this, but it doesn't hurt to have it
# anyway)
[[ $MANPATH =~ :$ ]] || export MANPATH=$MANPATH:
