######################################################################
#  This file is part of cddlib.
#
#        Copyright (C) 2020 Julian Rüth
#
#  cddlib is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  cddlib is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with cddlib. If not, see <https://www.gnu.org/licenses/>.
#####################################################################

import sys

try:
  input("Are you sure you are on the master branch which is identical to origin/master and the only pending changes are a version bump in the configure.ac of the library and possible author annotations in the news files and in AUTHORS? [ENTER]")
except KeyboardInterrupt:
  sys.exit(1)

sys.path.insert(0, 'recipe/snippets/rever')

import dist

$PROJECT = 'cddlib'

$ACTIVITIES = [
    'version_bump',
    'changelog',
    'dist',
    'tag',
    'push_tag',
    'ghrelease',
]

$VERSION_BUMP_PATTERNS = [
    ('configure.ac', r'AC_INIT', r'AC_INIT([cddlib], [$VERSION])'),
    ('lib-src/cddtypes.h', r'#define dd_DDVERSION', r'#define dd_DDVERSION   "Version $VERSION"'),
]

$CHANGELOG_FILENAME = 'ChangeLog'
$CHANGELOG_TEMPLATE = 'TEMPLATE.rst'
$PUSH_TAG_REMOTE = 'git@github.com:cddlib/cddlib.git'

$GITHUB_ORG = 'cddlib'
$GITHUB_REPO = 'cddlib'

$GHRELEASE_ASSETS = ['cddlib-' + $VERSION + '.tar.gz']
