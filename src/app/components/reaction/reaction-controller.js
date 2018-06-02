var host = "http://35.198.250.233/"
var RUN_REACTION_API = host + "upload?cleanInput=1"
var REACTION_API = host + "resource/reaction"
var GET_OUTPUT_API = host + 'download/'
var CHECK_SMART_API = host + 'tools/smart/'
var TYPE_REACTION = 0;
var TYPE_OPTION = 1;
class ReactionController {
    constructor($scope, $log, $http, $q, $timeout, $window, $mdDialog) {
        console.log("ReactionController")
        var self = this;
        self.simulateQuery = false;
        self.isDisabled = false;
        self.reagentRequired = false;
        self.selectedReactionList = [];
        self.selectedOptionList = [];
        self.loopChoice = createLoopChoice();
        // list of `state` value/display objects
        self.states = loadAll();
        self.querySearch = querySearch;
        self.selectedItemChange = selectedItemChange;
        self.searchTextChange = searchTextChange;
        self.newState = newState;
        self.sendRequest = sendRequest;
        self.addReaction = addReaction;
        self.addNewReaction = addNewReaction;
        self.editReaction = editReaction;
        self.removeReaction = removeReaction;
        self.removeOption = removeOption;
        self.downloadFile = downloadFile;
        self.outputReady = false;
        self.onDropFile = onDropFile;
        self.getCurrentWorkFlow = getCurrentWorkFlow;
        self.showDropZone = showDropZone;

        $("#sortable").sortable();
        $("#sortable").disableSelection();

        function editReaction(ev) {
            showDialog(ev, "EditReactionDialogController", 'app/components/reaction/dialog/edit-reaction-dialog.html')
        }
        function addNewReaction(ev) {
            showDialog(ev, "DialogController", 'app/components/reaction/dialog/new-reaction-dialog.html')
        }
        function showDialog(ev, controller, html) {
            $mdDialog.show({
                controller: controller + " as vm",
                templateUrl: html,
                parent: angular.element(document.body),
                targetEvent: ev,
                clickOutsideToClose: true,
                fullscreen: $scope.customFullscreen // Only for -xs, -sm breakpoints.
            })
        }
        function showDropZone(){
            return !(self.selectedOptionList.length >0 || self.selectedReactionList.length>0)
        }
        function addReaction(reaction) {
            if (reaction != null) {
                
                self.selectedReactionList.push(reaction);
            }
        }
        function addOption(option) {
            console.log("add option")
            console.log(option)
            if (option != null) {
                self.selectedOptionList.push(option);
            }
        }
        function removeReaction(id) {
            console.log("remove: " + id)
            for (var i = 0; i < self.selectedReactionList.length; i++) {
                if (self.selectedReactionList[i].id == id) {
                    self.selectedReactionList.splice(i, 1);
                    break;
                }
            }
        }
        function removeOption(id) {
            for (var i = 0; i < self.selectedOptionList.length; i++) {
                if (self.selectedOptionList[i].id == id) {
                    self.selectedOptionList.splice(i, 1);
                    break;
                }
            }
        }
        function downloadFile() {
            $window.open(GET_OUTPUT_API + self.outputFile, '_blank');
        }
      
        function getCurrentWorkFlow() {

            var reaction = [];
            var options = [];
            self.selectedReactionList.forEach(function (item) {
                reaction.push(item)
            });
            self.selectedOptionList.forEach(function (item) {
                options.push(item);
            });
            var workFlow = { reaction, options };
            var json = JSON.stringify(workFlow);
            var properties = { type: 'application/json' };
            var file;
            var data = [];
            data.push(json)
            var filename = "config.json";
            try {
                file = new File(data, filename, properties);
            } catch (e) {
                file = new Blob(data, properties);
            }
            var url = window.URL.createObjectURL(file);
            var a = document.getElementById("config")
            a.href = url;
            a.download = filename;
            a.click();
            window.URL.revokeObjectURL(url);
            return JSON.stringify(workFlow);
        }
        function sendRequest() {
            var payload = new FormData();
            payload.append("file", $scope.file);
            var reaction = [];
            var options = [];
            self.selectedReactionList.forEach(function (item) {
                reaction.push({ id: item.id, reagent: item.input_reagent, loop: item.loop })
            });
            self.selectedOptionList.forEach(function (item) {
                let option = Object.assign({}, item);
                option.template = "";
                option.warn = "";
                options.push(option); 
            });
            payload.append("reaction", JSON.stringify(reaction));
            payload.append("options", JSON.stringify(options));
            self.sampleResult = "";
            $http({
                url: RUN_REACTION_API,
                method: 'POST',
                data: payload,
                //assign content-type as undefined, the browser
                //will assign the correct boundary for us
                headers: { 'Content-Type': undefined },
                //prevents serializing payload.  don't do it.
                transformRequest: angular.identity
            }).then(function (response) {
                console.log(response);
                try {
                    var data = response.data;
                    console.log(data)
                    self.outputFile = data.output;
                    self.outputReady = true;
                    if (typeof data.sample != 'undefined') {
                        self.sampleResult = data.sample;
                    }
                } catch (error) {
                    console.log(error);
                }

            }, function (error) {
                console.log(error)
            })

        }
        function loadAll() {
            var option_dir = "app/components/reaction/options/";
            return new Promise(function (resolve, reject) {
                $http.get(REACTION_API).then(function (response) {
                    self.reactionList = response.data.reaction;
                    self.reactionList.forEach(function (element) {
                        element.display = element.name + " " + element.smart;
                        element.loop = 1;
                        element.type = TYPE_REACTION;
                    });
                    self.reactionList.push({ id: 0, display: "3D Maker", type: TYPE_OPTION, ionize:false, pH:0, description: "", removeSalt:true, warn: "This option will take time to process", template:option_dir+"3dmaker.html" })
                    self.reactionList.push({ id: 1, display: "Preparator", type: TYPE_OPTION, description: "", ionize:false, pH:0, addHs:true, removeSalt:true, template:option_dir+"preparator.html" })
                    resolve()
                });
            })
        }
        function querySearch(query) {
            //console.log(self.reactionList)
            if (typeof self.reactionList == 'undefined') {
                var deferred = $q.defer();
                loadAll().then((success) => {
                    var results = query ? self.reactionList.filter(createFilterFor(query)) : self.reactionList;
                    deferred.resolve(results)
                })
                return deferred.promise;
            } else {
                var results = query ? self.reactionList.filter(createFilterFor(query)) : self.reactionList;
                return results
            }

        }
        function searchTextChange(text) {
            $log.info('Text changed to ' + text);
        }

        function selectedItemChange(item) {
            if (item != null) {
                switch (item.type) {
                    case TYPE_REACTION:
                        addReaction(item);
                        break;
                    case TYPE_OPTION:
                        addOption(item);
                        break;
                }
                $log.info('Item changed to ' + JSON.stringify(item));
            }
        }
        function newState(state) {
            alert("Sorry! You'll need to create a Constitution for " + state + " first!");
        }
        function createFilterFor(query) {
            console.log("create filter");
            var lowercaseQuery = query.toLowerCase();

            return function filterFn(state) {
                return (state.display.toLowerCase().includes(lowercaseQuery));
            };
        }
        function createLoopChoice() {
            var choice = [];
            for (var i = 1; i <= 5; i++) {
                choice.push({ display: i, value: i });
            }
            return choice;
        }
       
        function findBlock(id, type){
            for(var i = 0;i<self.reactionList.length;i++){
                let item = self.reactionList[i];
                if(item.id==id && item.type == type){
                    return item;
                }
            }
            return null;
        }
        function onDropFile(data) {
            var object = JSON.parse(data);
            var reaction = object.reaction;
            console.log(object)
            reaction.forEach(function(item){
                self.addReaction(item);
            })
           
           
            var option_list = object.options;
            option_list.forEach(function(item){
                addOption(item)
            });
            $scope.$apply();
          
        }


    }

}

angular
    .module('app')
    .component('reaction', {
        templateUrl: 'app/components/reaction/index.html',
        controller: ReactionController,
        controllerAs: "vm"
    })
    .directive('chooseFile', function () {
        return {
            link: function (scope, elem, attrs) {
                var button = elem.find('button');
                var input = angular.element(elem[0].querySelector('input#fileInput'));
                button.bind('click', function () {
                    input[0].click();
                });
                input.bind('change', function (e) {
                    scope.$apply(function () {
                        var files = e.target.files;
                        if (files[0]) {
                            scope.file = files[0];
                            scope.fileName = files[0].name;

                        } else {
                            scope.fileName = null;
                        }
                    });
                });
            }
        };
    })

