#!/bin/tcsh

foreach x ( inputs/*.yaml )
    grep -c $x *.yaml
    if ( $? != 0 ) then
	echo $x
    endif
end
