#!/bin/bash
touch $2
for jd in `cat $1`
do
  curl --referer http://astroutils.astronomy.ohio-state.edu/time/bjd2utc.html -d jds=$jd \
  -d ra=14.28275000 \
  -d raunits=hours -d dec=-49.94500000 \
  -d spaceobs=none "http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.php" \
   2>/dev/null | sed -n '/^[0-9]/ p' >> $2
done
