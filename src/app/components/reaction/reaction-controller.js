var host = "http://35.198.250.233/"
var RUN_REACTION_API = host + "upload"
var GET_REACTION_API = host + "resource/reaction"
class ReactionController {
    constructor($scope, $log, $http, $q, $timeout) {
        console.log("ReactionController")
        var self = this;
        self.simulateQuery = false;
        self.isDisabled = false;
        self.message = "hello"
        self.reagentRequired = false;
        self.selectedReactionList = [];
        // list of `state` value/display objects
        // self.states = loadAll();
        self.querySearch = querySearch;
        self.selectedItemChange = selectedItemChange;
        self.searchTextChange = searchTextChange;
        self.sendRequest = sendRequest;
        self.newState = newState;
        self.addReaction = addReaction;
        self.removeReaction = removeReaction;
        self.currentSearchedReaction = {};
        $("#sortable").sortable();
        $("#sortable").disableSelection();
        function checkReagentInput() {
            self.reagentRequired = false;
            self.selectedReactionList.forEach(function (item) {
                if (item.reagent) {
                    self.reagentRequired = true;
                }
            });
        }
        function addReaction() {
            self.selectedReactionList.push(self.currentSearchedReaction);
            checkReagentInput();
        }
        function removeReaction(id) {
            console.log("remove: " + id)
            for (var i = 0; i < self.selectedReactionList.length; i++) {
                if (self.selectedReactionList[i].id == id) {
                    self.selectedReactionList.splice(i, 1);
                    break;
                }
            }
            checkReagentInput();
        }
        function sendRequest() {

            var payload = new FormData();
            payload.append("file", $scope.file);
            var reaction = [];
            self.selectedReactionList.forEach(function (item) {
                reaction.push({ reaction: item.id, reagent: item.input_reagent })

            });
            payload.append("reaction", JSON.stringify(reaction));

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
            }, function (error) {
                console.log(error)
            })

        }
        function loadAll() {
            return new Promise(function (resolve, reject) {
                $http.get(GET_REACTION_API).then(function (response) {
                    self.reactionList = response.data.reaction;
                    self.reactionList.forEach(function (element) {
                        element.display = element.name + " " + element.smart;
                    });
                    resolve()
                });
            })
        }
        function querySearch(query) {
            //console.log(self.reactionList)
            if (typeof self.reactionList == 'undefined') {
                console.log("come in here")


                var deferred = $q.defer();
                console.log("create differ");
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
            self.currentSearchedReaction = item;
            $log.info('Item changed to ' + JSON.stringify(item));

        }
        function newState(state) {
            alert("Sorry! You'll need to create a Constitution for " + state + " first!");
        }
        function createFilterFor(query) {
            console.log("create filter");
            var lowercaseQuery = query.toLowerCase();

            return function filterFn(state) {

                return (state.name.includes(lowercaseQuery));
            };

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
    });