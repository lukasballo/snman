Turn costs in R5:
- public static final int LEFT_TURN = 30;
- public static final int STRAIGHT_ON = 0;
- public static final int RIGHT_TURN = 10;
- public static final int U_TURN = 90; // penalize U turns extremely heavily

Turn costs for various types of turns, in seconds. These are all specified for drive-on-right countries,
and reversed in drive-on-left countries

    public int computeTurnCost (int fromEdge, int toEdge, StreetMode streetMode) {
        if (streetMode == StreetMode.CAR) {
            double angle = calculateNewTurnAngle(fromEdge, toEdge);

            if (angle < 27)
                return STRAIGHT_ON;
            else if (angle < 153)
                return driveOnRight ? LEFT_TURN : RIGHT_TURN;
            else if (angle < 207)
                return U_TURN;
            else if (angle < 333)
                return driveOnRight ? RIGHT_TURN : LEFT_TURN;
            else
                return STRAIGHT_ON;
        }

        return 0;
    }

https://github.com/mcheung610/r5/blob/master/src/main/java/com/conveyal/r5/streets/TurnCostCalculator.java

