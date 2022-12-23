# chembalancer
A common task in chemistry class is to balance a chemical reaction.Typically, the substances on each side of the reaction are given and the challenge is to find integer coefficients in front of every substance to make the reaction balanced. This code does just that!!!

## Example usecase:
```console
$pipenv run python ./chembalancer.py
```
Output:
```
Enter reaction formula:
C6H12O6 + O2 = H2O + CO2
Reaction is uniquely balanced
C6H12O6 + 6 O2 -> 6 H2O + 6 CO2
```

## Setup
1. Install `pipenv`
2. Install dependencies by running `pipenv install`
3. Run the chembalancer script. For example by running  `pipenv run python ./chembalancer.py`

## Details
The script will parse input according to the following assumptions. The substances are listed using notations from the periodic table, the substances are separated by '+'-signs and the two sides of the reaction is separated by an '='-sign.

Thereafter the problem is converted into a system of linear integer equations, which is solved by bringing it to [Smith Normal Form](https://en.wikipedia.org/wiki/Smith_normal_form). In general the solution to this system may not be unique. In the case of several linearly independant solutions an integer approach to the [Conflict Resolution](https://www.semanticscholar.org/paper/Conflict-Resolution-Korovin-Tsiskaridze/c1b16de4d26b6efe97d195e7a85ac36377badba2) algorithm is used to solve the set of inequalities required to find solutions with only positive coefficients.
