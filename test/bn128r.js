const fieldTests = require("./tester/fieldtests.js");

const bigInt = require("big-integer");

fieldTests(new bigInt("21888242871839275222246405745257275088548364400416034343698204186575808495617"), "bn128r");
