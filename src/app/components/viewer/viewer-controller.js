class ViewerController {
  constructor($scope, $http) {
    console.log("ViewerController")
    var self = this;
    self.reload = reload;
    self.executeFile = executeFile;
    self.sdf_list = [];
    self.showIndex = 0;
    self.molecule_name = "";
    self.changeIndex = changeIndex;
    self.hasOutput = false;

    var reader = new FileReader();
    reader.onload = function (event) {
      var data = event.target.result;
      self.sdf_list = data.split("$$$$")
      self.sdf_list = self.sdf_list.map(s => s.trim());

      //console.log(self.sdf_list[1])
      self.hasOutput = true;
      changeIndex(0);
      $scope.$apply();

    };
    $scope.$on('update_sdf_sample', function (event, args) {
      console.log("update sdf form viewController");
      console.log(args);
      var sdf_file = args;
      var url = GET_RAW_SDF + sdf_file;
      $http({
        url: url,
        method: 'GET',
        headers: {
          'Content-Type': "text/plain"
        },
        transformRequest: angular.identity
      }).then(function (response) {
        var data = response.data;
        self.sdf_list = data.split("$$$$")
        self.sdf_list = self.sdf_list.map(s => s.trim());
        changeIndex(0);

      }, function (error) {
        console.log(error)
      })

      // do what you want to do
    });
    function executeFile() {
      console.log($scope.file)
      reader.readAsText($scope.file);
    }

    var glmol01 = new GLmol('glmol01', true);



    function loadFile() {
      var file = $('#glmol01_file').get(0);
      if (file) file = file.files;
      if (!file || !window.FileReader || !file[0]) {
        alert("No file is selected. Or File API is not supported in your browser. Please try Firefox or Chrome.");
        return;
      }
      $('#loading').show();
      var reader = new FileReader();
      reader.onload = function () {
        $('#glmol01_src').val(reader.result);
        glmol01.loadMolecule();
        $('#loading').hide();
      };
      reader.readAsText(file[0]);
    }

    function saveImage() {
      glmol01.show();
      var imageURI = glmol01.renderer.domElement.toDataURL("image/png");
      window.open(imageURI);
    }
    function changeIndex(value) {
      self.showIndex += value;
      if (self.showIndex <= 0) {
        self.showIndex += self.sdf_list.length -1;
      }
      self.showIndex %= self.sdf_list.length;
      console.log("change index " + self.showIndex);

      $("#glmol01_src").val(self.sdf_list[self.showIndex]);
      // console.log(self.sdf_list[self.showIndex]);
      reload();
    }
    function reload() {

      glmol01.loadMolecule();
      console.log("reload: " + glmol01.molecule_name)
      self.molecule_name = glmol01.molecule_name;
      // $scope.$apply();
      // console.log( $("#glmol01_src").val().split("$$$$"))
      glmol01.defineRepresentation = defineRepFromController;
      glmol01.rebuildScene();
      glmol01.show();
    }


    function defineRepFromController() {
      var idHeader = "#" + this.id + '_';

      var time = new Date();
      var all = this.getAllAtoms();
      if ($(idHeader + 'biomt').attr('checked') && this.protein.biomtChains != "") all = this.getChain(all, this.protein.biomtChains);
      var allHet = this.getHetatms(all);
      var hetatm = this.removeSolvents(allHet);

      console.log("selection " + (+new Date() - time));
      time = new Date();

      this.colorByAtom(all, {});
      var colorMode = "chainbow"
      if (colorMode == 'ss') {
        this.colorByStructure(all, 0xcc00cc, 0x00cccc);
      } else if (colorMode == 'chain') {
        this.colorByChain(all);
      } else if (colorMode == 'chainbow') {
        this.colorChainbow(all);
      } else if (colorMode == 'b') {
        this.colorByBFactor(all);
      } else if (colorMode == 'polarity') {
        this.colorByPolarity(all, 0xcc0000, 0xcccccc);
      }
      console.log("color " + (+new Date() - time));
      time = new Date();

      var asu = new THREE.Object3D();
      var mainchainMode = "thickRibbon" //mod
      var showMainChain = true; //mod
      var doNotSmoothen = false //mod
      if (showMainChain) {
        if (mainchainMode == 'ribbon') {
          this.drawCartoon(asu, all, doNotSmoothen);
          this.drawCartoonNucleicAcid(asu, all);
        } else if (mainchainMode == 'thickRibbon') {
          this.drawCartoon(asu, all, doNotSmoothen, this.thickness);
          this.drawCartoonNucleicAcid(asu, all, null, this.thickness);
        } else if (mainchainMode == 'strand') {
          this.drawStrand(asu, all, null, null, null, null, null, doNotSmoothen);
          this.drawStrandNucleicAcid(asu, all);
        } else if (mainchainMode == 'chain') {
          this.drawMainchainCurve(asu, all, this.curveWidth, 'CA', 1);
          this.drawMainchainCurve(asu, all, this.curveWidth, 'O3\'', 1);
        } else if (mainchainMode == 'cylinderHelix') {
          this.drawHelixAsCylinder(asu, all, 1.6);
          this.drawCartoonNucleicAcid(asu, all);
        } else if (mainchainMode == 'tube') {
          this.drawMainchainTube(asu, all, 'CA');
          this.drawMainchainTube(asu, all, 'O3\''); // FIXME: 5' end problem!
        } else if (mainchainMode == 'bonds') {
          this.drawBondsAsLine(asu, all, this.lineWidth);
        }
      }
      var line = false; //mod
      if (line) {
        this.drawBondsAsLine(this.modelGroup, this.getSidechains(all), this.lineWidth);
      }
      console.log("mainchain " + (+new Date() - time));
      time = new Date();
      var showBased = true;
      if (showBased) {
        var hetatmMode = "nuclStick"
        if (hetatmMode == 'nuclStick') {
          this.drawNucleicAcidStick(this.modelGroup, all);
        } else if (hetatmMode == 'nuclLine') {
          this.drawNucleicAcidLine(this.modelGroup, all);
        } else if (hetatmMode == 'nuclPolygon') {
          this.drawNucleicAcidLadder(this.modelGroup, all);
        }
      }

      var target = $(idHeader + 'symopHetatms').attr('checked') ? asu : this.modelGroup;
      var showNonBond = false;
      if (showNonBond) {
        var nonBonded = this.getNonbonded(allHet);
        var nbMode = "nb_sphere";
        if (nbMode == 'nb_sphere') {
          this.drawAtomsAsIcosahedron(target, nonBonded, 0.3, true);
        } else if (nbMode == 'nb_cross') {
          this.drawAsCross(target, nonBonded, 0.3, true);

        }
      }
      var showhetatm = true;
      if (showhetatm) {
        var hetatmMode = "ballAndStick2"
        if (hetatmMode == 'stick') {
          this.drawBondsAsStick(target, hetatm, this.cylinderRadius, this.cylinderRadius, true);
        } else if (hetatmMode == 'sphere') {
          this.drawAtomsAsSphere(target, hetatm, this.sphereRadius);
        } else if (hetatmMode == 'line') {
          this.drawBondsAsLine(target, hetatm, this.curveWidth);
        } else if (hetatmMode == 'icosahedron') {
          this.drawAtomsAsIcosahedron(target, hetatm, this.sphereRadius);
        } else if (hetatmMode == 'ballAndStick') {
          this.drawBondsAsStick(target, hetatm, this.cylinderRadius / 2.0, this.cylinderRadius, true, false, 0.3);
        } else if (hetatmMode == 'ballAndStick2') {
          this.drawBondsAsStick(target, hetatm, this.cylinderRadius / 2.0, this.cylinderRadius, true, true, 0.3);
        }

      }
      console.log("hetatms " + (+new Date() - time));
      time = new Date();

      var projectionMode = "perspective" //mods
      if (projectionMode == 'perspective') this.camera = this.perspectiveCamera;
      else if (projectionMode == 'orthoscopic') this.camera = this.orthoscopicCamera;
      var bgcolor = "0x888888"
      this.setBackground(parseInt(0));
      var cell = false; //mod
      if (cell) {
        this.drawUnitcell(this.modelGroup);
      }
      var biomt = false; //mod
      if (biomt) {
        this.drawSymmetryMates2(this.modelGroup, asu, this.protein.biomtMatrices);
      }
      var packing = false; //mod
      if (packing) {
        this.drawSymmetryMatesWithTranslation2(this.modelGroup, asu, this.protein.symMat);
      }
      this.modelGroup.add(asu);
    };

    glmol01.defineRepresentation = defineRepFromController;

  }
}
angular
  .module('app')
  .controller('ViewerController', ViewerController)
  .component('viewer', {
    templateUrl: 'app/components/viewer/index.html',
    controller: ViewerController,
    controllerAs: "vm"
  })
