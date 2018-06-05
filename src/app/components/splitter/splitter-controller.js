class SplitterController {
    constructor($scope) {
        var self = this;
        self.data = null;
        self.executeFile = executeFile;
        self.splitList = [];
        self.fileList = [];
        self.filename = "output";
        self.itemPerFile = 1000;
        self.autoDownload = false;
        self.number_choice = [100, 300, 500, 800, 1000];
        var EXTENSION = { smile: { name: 'smi', delimiter: '\n' }, sdf: { name: 'sdf', delimiter: '$$$$' } };
        console.log("SplitterController")
        var reader = new FileReader();

        reader.onload = function (event) {
            self.data = event.target.result;
            //console.log(self.data)
            splitFile(self.itemPerFile);
            $scope.$apply();
        };
        function executeFile() {
            console.log($scope.file);
            self.filename = $scope.file.name;
            var parts = self.filename.split(".");
            self.no_extension_filename = parts[0];
            self.extension = parts[1];
            reader.readAsText($scope.file);
        }
        function splitFile(item_number) {
            if (self.data == null || typeof self.data == 'undefined') return;
            var properties = {
                type: 'text/plain'
            };
            var itemList = [];
            var extension = EXTENSION.smile;
            if (self.extension == EXTENSION.smile.name) {
                extension = EXTENSION.smile;
            } else if (self.extension == EXTENSION.sdf.name) {
                extension = EXTENSION.sdf;
            }
            itemList = self.data.split(extension.delimiter);
            self.splitList = [];
            var str = "";
            for (var i = 0; i < itemList.length; i++) {
                str += itemList[i] + extension.delimiter;
                if (i % item_number == 1 || i == itemList.length - 1) {
                    self.splitList.push(str);
                    str = "";
                }
            }
            self.fileList = [];
            for (var i = 0; i < self.splitList.length; i++) {
                var file;
                var filename = self.no_extension_filename + "_" + (i+1) + "." + self.extension;
                try {
                    file = new File([self.splitList[i]], filename, properties);
                } catch (e) {
                    file = new Blob([self.splitList[i]], properties);
                }
                var url = window.URL.createObjectURL(file);
                self.fileList.push({ url, filename });
            }
            $scope.$apply();

            var a_list = document.getElementsByClassName("output_files")
            console.log(a_list);
            for (var i = 0; i < a_list.length; i++) {
                a_list[i].href = self.fileList[i].url;
                a_list[i].download = self.fileList[i].filename
                if (self.autoDownload) {
                    a_list[i].click();
                }
                window.URL.revokeObjectURL(self.fileList[i].url);
            }

            return true;
        }
    }

}
angular
    .module('app')
    .controller('SplitterController', SplitterController)
    .component('splitter', {
        templateUrl: 'app/components/splitter/index.html',
        controller: SplitterController,
        controllerAs: "vm"
    })
