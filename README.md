#### Build

1. Install dependencies

```bash
git clone https://github.com/jplumail/TiPi
cd TiPi
mvn install
```

2. Build JAR

Clone this repo, and do:
```bash
mvn package
```

The JAR file `target/microTiPi-1.0.jar` is executable:
```bash
java -jar target/microTiPi-1.0.jar -help
```

#### Usage

There are 2 subcommands:
- deconv
- blinddeconv

See the help of the commands:

```bash
java -jar microTiPi-1.0.jar deconv -help
java -jar microTiPi-1.0.jar blinddeconv -help
```

The `deconv` command outputs a deconvolved image.

The `blinddeconv` command outputs a deconvolved image and the psf found.