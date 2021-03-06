#!/usr/bin/env node

const bigInt=require("big-integer");
const path = require("path");
const util = require("util");
const renderFile = util.promisify(require("ejs").renderFile);

const runningAsScript = !module.parent;

const montgomeryBuilder = require("./montgomerybuilder");

class ZqBuilder {
    constructor(q, name) {
        const self = this;
        this.q=bigInt(q);
        this.n64 = Math.floor((this.q.bitLength() - 1) / 64)+1;
        this.name = name;
        this.bigInt = bigInt;
        this.lastTmp=0;
        this.global = {};
        this.global.tmpLabel = function(label) {
            self.lastTmp++;
            label = label || "tmp";
            return label+"_"+self.lastTmp;
        };
        this.montgomeryBuilder = montgomeryBuilder;
    }

    constantElement(v) {
        let S = "";
        const mask = bigInt("FFFFFFFFFFFFFFFF", 16);
        for (let i=0; i<this.n64; i++) {
            if (i>0) S = S+",";
            let shex = v.shiftRight(i*64).and(mask).toString(16);
            while (shex.length <16) shex = "0" + shex;
            S = S + "0x" + shex;
        }
        return S;
    }

}

async function buildField(q, name) {
    const builder = new ZqBuilder(q, name);

    let asm = await renderFile(path.join(__dirname, "fr.asm.ejs"), builder);
    const cpp = await renderFile(path.join(__dirname, "fr.cpp.ejs"), builder);
    const hpp = await renderFile(path.join(__dirname, "fr.hpp.ejs"), builder);

    return {asm: asm, hpp: hpp, cpp: cpp};
}

if (runningAsScript) {
    const fs = require("fs");
    var argv = require("yargs")
        .usage("Usage: $0 -q [primeNum] -n [name] -oc [out .c file] -oh [out .h file] -oa [out .asm file]")
        .demandOption(["q","n"])
        .alias("q", "prime")
        .alias("n", "name")
        .argv;

    const q = bigInt(argv.q);

    const asmFileName =  (argv.oa) ? argv.oa : argv.name.toLowerCase() + ".asm";
    const hFileName =  (argv.oh) ? argv.oh : argv.name.toLowerCase() + ".hpp";
    const cFileName =  (argv.oc) ? argv.oc : argv.name.toLowerCase() + ".cpp";


    buildField(q, argv.name).then( (res) => {
        fs.writeFileSync(asmFileName, res.asm, "utf8");
        fs.writeFileSync(hFileName, res.hpp, "utf8");
        fs.writeFileSync(cFileName, res.cpp, "utf8");
    });

} else {
    module.exports = buildField;
}
