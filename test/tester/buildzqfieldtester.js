const chai = require("chai");
const assert = chai.assert;

const fs = require("fs");
var tmp = require("tmp-promise");
const path = require("path");
const util = require("util");
const exec = util.promisify(require("child_process").exec);

const buildZqField = require("../../index.js").buildZqField;

module.exports = testField;

// Compiled testers, one per (prime, mode). The binary does not depend on the test vector.
const testers = new Map();

async function buildTester(prime, mode) {
    const dir = await tmp.dir({prefix: "ffiasm_", unsafeCleanup: true });

    const source = await buildZqField(prime, "Fr");

    // console.log(dir.path);

    const files = {
        "fr.asm": source.asm,
        "fr.hpp": source.hpp,
        "fr.cpp": source.cpp,
        "fr_element.hpp": source.element_hpp,
        "fr_generic.cpp": source.generic_cpp,
        "fr_raw_generic.cpp": source.raw_generic_cpp,
        "fr_raw_arm64.s": source.raw_arm64_s
    };
    for (const name of Object.keys(files)) {
        await fs.promises.writeFile(path.join(dir.path, name), files[name], "utf8");
    }

    await exec(`cp  ${path.join(__dirname,  "tester.cpp")} ${dir.path}`);

    let flags = "";
    let sources;
    if (mode == "generic") {
        sources = ["tester.cpp", "fr.cpp", "fr_generic.cpp", "fr_raw_generic.cpp"];
    } else if (process.arch == "x64") {
        if (process.platform === "darwin") {
            await exec("nasm -fmacho64 --prefix _ " +
                ` ${path.join(dir.path,  "fr.asm")}`
            );
        }  else if (process.platform === "linux") {
            await exec("nasm -felf64 " +
                ` ${path.join(dir.path,  "fr.asm")}`
            );
        } else throw("Unsupported platform");
        flags = " -DUSE_ASM -DARCH_X86_64";
        sources = ["tester.cpp", "fr.cpp", "fr.o"];
    } else if (process.arch == "arm64") {
        flags = " -DUSE_ASM -DARCH_ARM64";
        sources = ["tester.cpp", "fr.cpp", "fr_generic.cpp", "fr_raw_generic.cpp", "fr_raw_arm64.s"];
    } else throw("Unsupported architecture");

    // Homebrew on Apple Silicon installs gmp out of the default search paths
    if ((process.platform === "darwin")&&(fs.existsSync("/opt/homebrew/include/gmp.h"))) {
        flags += " -I/opt/homebrew/include -L/opt/homebrew/lib";
    }

    await exec("g++" + flags +
               sources.map( (f) => " " + path.join(dir.path, f)).join("") +
               ` -o ${path.join(dir.path, "tester")}` +
               " -lgmp -g"
    );

    return dir;
}

async function  testField(prime, test, mode) {
    tmp.setGracefulCleanup();

    mode = mode || "asm";
    const key = prime.toString() + "_" + mode;
    if (!testers.has(key)) testers.set(key, buildTester(prime, mode));
    const dir = await testers.get(key);

    const inLines = [];
    for (let i=0; i<test.length; i++) {
        for (let j=0; j<test[i][0].length; j++) {
            inLines.push(test[i][0][j]);
        }
    }
    inLines.push("");

    await fs.promises.writeFile(path.join(dir.path, "in.tst"), inLines.join("\n"), "utf8");

    await exec(`${path.join(dir.path, "tester")}` +
        ` <${path.join(dir.path, "in.tst")}` +
        ` >${path.join(dir.path, "out.tst")}`);

    const res = await fs.promises.readFile(path.join(dir.path, "out.tst"), "utf8");
    const resLines = res.split("\n");

    for (let i=0; i<test.length; i++) {
        const expected = test[i][1].toString();
        const calculated = resLines[i];

        if (calculated != expected) {
            console.log("FAILED");
            for (let j=0; j<test[i][0].length; j++) {
                console.log(test[i][0][j]);
            }
            console.log("Should Return: " + expected);
            console.log("But Returns: " + calculated);
        }

        assert.equal(calculated, expected);
    }

}

